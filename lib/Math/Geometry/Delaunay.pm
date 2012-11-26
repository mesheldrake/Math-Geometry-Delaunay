package Math::Geometry::Delaunay;

use 5.008;
use warnings;
use strict;
use Carp qw(carp);;
use Exporter();

our @ISA = qw(Exporter);
our $VERSION;

BEGIN {
    use XSLoader;
    $VERSION = '0.04';
    XSLoader::load('Math::Geometry::Delaunay');
    }

use constant {
    TRI_CONSTRAINED  => 'Y',
    TRI_CONFORMING   => 'Dq0',
    TRI_CCDT         => 'q',
    TRI_VORONOI      => 'v'
    };

our @EXPORT_OK = qw(TRI_CONSTRAINED TRI_CONFORMING TRI_CCDT TRI_VORONOI);
our @EXPORT = qw();

sub new {
    my $class = shift;
    my $self = {};
    $self->{in}     = Math::Geometry::Delaunay::Triangulateio->new();
    $self->{out}    = undef;
    $self->{vorout} = undef;
    $self->{poly}   = {
	regions      => [],
	holes        => [],
	polylines    => [],
	points       => [],
	segments     => [],
        outnodes     => [], #for cache, first time C output arrays are imported
        voutnodes    => [], #for cache
        segptrefs    => {},  #used to avoid dups of points added with addSegment
        };

    $self->{optstr} = '';
    # Triangle switches
    # prq__a__uAcDjevngBPNEIOXzo_YS__iFlsCQVh
    # where __ is an optional number
    $self->{a} = -1; # max tri area
    $self->{q} = -1; # quality min angle
    $self->{e} = 0; # produce edges switch
    $self->{v} = 0; # voronoi switch
    $self->{n} = 0; # neighbors switch
    $self->{N} = 0; # suppress node output
    $self->{E} = 0; # suppress element output
    $self->{O} = 0; # suppress holes - ignore input holes
    $self->{o2}= 0; # subparametric switch (for 6 pts/tri)
    $self->{Q} = 1; # quiet - switch for Triangle's printed output    
    $self->{V} = 0; # verbose - from 0 to 3 for increasing verbosity

    bless $self, $class;
    return $self;
    }

# triangulatio interfaces
sub in     {return $_[0]->{in};}
sub out    {return $_[0]->{out};}
sub vorout {return $_[0]->{vorout};}

# getter/setters for the triangulate switches that take numbers

sub area_constraint { # -1 for disabled
    if (@_>1) {$_[0]->{a}=$_[1];}
    return $_[0]->{a};
    }
sub minimum_angle { # "q" switch, in degrees, -1 for disabled
    if (@_>1) {$_[0]->{q}=$_[1];}
    return $_[0]->{q};
    }
sub subparametric {
    if (@_>1) {$_[0]->{o2}=$_[1]?1:0;}
    return $_[0]->{o2};
    }
sub doEdges {
    if (@_>1) {$_[0]->{e}=$_[1]?1:0;}
    return $_[0]->{e};
    }
sub doVoronoi {
    if (@_>1) {$_[0]->{v}=$_[1]?1:0;}
    return $_[0]->{v};
    }
sub doNeighbors {
    if (@_>1) {$_[0]->{n}=$_[1]?1:0;}
    return $_[0]->{n};
    }
sub quiet {
    if (@_>1) {$_[0]->{Q}=$_[1]?1:0;}
    return $_[0]->{Q};
    }
sub verbose { # 0 to 3
    if (@_>1) {$_[0]->{V}=$_[1]?1:0;}
    return $_[0]->{V};
    }

# everything to add input geometry

sub addRegion {
    my $self = shift;
    my $poly = shift;
    my $attribute = @_ ? shift : undef;
    my $area = @_ ? shift:undef;
    my $point_inside = @_ ? shift : undef; # not expected, but we'll use it

    if (@{$poly}==1) {
        carp "first arg to addRegion should be a polygon, or point";
        return;
        }
    elsif (@{$poly}==2 && !ref($poly->[0])) { # a region identifying point
        $point_inside = $poly;
        }
    else {
        $self->addPolygon($poly);
        }

    my $ray; # return ray used for $point_inside calc for debugging, for now

    if (!$point_inside) {
        ($point_inside, $ray) = get_point_in_polygon($poly);
        }
    if (defined $point_inside) {
        push @{$self->{poly}->{regions}}, [ $point_inside, $attribute, ($area && $area > 0) ? $area : -1 ];
        }
    return $point_inside, $ray;
    }

sub addHole {
    my $self = shift;
    my $poly = shift;
    my $point_inside = @_ ? shift : undef; # not expected, but we'll use it if available

    if (@{$poly}==1) {
        carp "first arg to addHole should be a polygon, or point";
        return;
        }
    elsif (@{$poly}==2 && !ref($poly->[0])) { # it's really the hole identifying point
        $point_inside = $poly;
        }
    else {
        $self->addPolygon($poly);
        }

    my $ray; # return ray used for $point_inside calc for debugging, for now

    if (!$point_inside) {
        ($point_inside, $ray) = get_point_in_polygon($poly);
        }
    if (defined $point_inside) {
        push @{$self->{poly}->{holes}}, $point_inside;
        }
    return $point_inside, $ray;
    }

sub addPolygon {
    my $self = shift;
    my $poly = shift;
    if    (@{$poly} == 1 ) {return $self->addPoints([$poly->[0]]);}
    push @{$self->{poly}->{polylines}}, ['polygon',$poly];
    return;
    }

sub addPolyline {
    my $self = shift;
    my $poly = shift;
    if    (@{$poly} == 1 ) {return $self->addPoints([$poly->[0]]);}
    push @{$self->{poly}->{polylines}}, ['polyline',$poly];
    return;
    }

sub addSegment {
    my $self = shift;
    my $points = shift;
    push @{$self->{poly}->{points}}, @{$points};
    return;
    }

sub addPoints { # points unaffiliated with PLSG segments
    my $self = shift;
    my $points = shift;
    push @{$self->{poly}->{points}}, @{$points};
    return;
    }

# compile all the input geometry in to Triangle-format lists
# set up option strings
# and initialize output lists

sub prepPoly {
    my $self = shift;
    my $optstr = shift;
    $self->{optstr} = '';
    # option string options:
    # prq__a__uAcDjevngBPNEIOXzo_YS__iFlsCQVh
    # where __ is an optional number
    $self->{optstr} .= ''.
                       $optstr.
                       ($self->{q} > -1?'q'.$self->{q}:'').
                       ($self->{a} > -1?'a'.$self->{a}:'').
                       ($self->{e}     ?'e':'').
                       ($self->{v}     ?'v':'').
                       ($self->{n}     ?'n':'').
                       'z'. # always number everything starting with zero
                       ($self->{o2}    ?'o2':'').
                       ($self->{Q}     ?'Q':'').
                       ($self->{V} >  0?(($self->{V} > 2) ? 'VVV' : ($self->{V} x 'V')) : '')
                       ;
    my @allpts;
    my @allsegs;

    if (@{$self->{poly}->{segments}}) {
        # The list of segments is the most likely to have duplicate points.
        # Some algorithms in this space result in lists of segments,
        # perhaps listing subsegments of intersecting segments,
        # or representing a boundary or polygon with out-of-order,
        # non-contiguous segment lists, where shared vertices are
        # repeated in each segment's record.
        # addSegments() is meant for that kind of data
        # and this is where we boil the duplicate points down,
        # so Triangle doesn't have to. 
        foreach my $seg (@{$self->{poly}->{segments}}) {
            if (!defined $self->{segrefs}->{$seg->[0]} &&
                !defined $self->{segrefs}->{$seg->[0]->[0]}.','.$self->{segrefs}->{$seg->[0]->[1]}
                ) {
                push @allpts, $seg->[0];
                $self->{segrefs}->{$seg->[0]}=$#allpts;
                $self->{segrefs}->{$seg->[0]->[0]}.','.$self->{segrefs}->{$seg->[0]->[1]} = $#allpts;
                }
            push @allsegs, $self->{segrefs}->{$seg->[0]};
            if (!defined $self->{segrefs}->{$seg->[1]} &&
                !defined $self->{segrefs}->{$seg->[1]->[0]}.','.$self->{segrefs}->{$seg->[1]->[1]}
                ) {
                push @allpts, $seg->[1];
                $self->{segrefs}->{$seg->[1]}=$#allpts;
                $self->{segrefs}->{$seg->[1]->[0]}.','.$self->{segrefs}->{$seg->[1]->[1]} = $#allpts;
                }
            push @allsegs, $self->{segrefs}->{$seg->[1]};            
            }
        }
    undef $self->{segrefs};

    if (@{$self->{poly}->{polylines}}) {
        # doing PSLG - add poly flag to options
        $self->{optstr} = 'p'.$self->{optstr};
        #set up points and segments lists for each polygon and polyline
        foreach my $polycont (@{$self->{poly}->{polylines}}) {
            my $poly = $polycont->[1];
            push @allpts, $poly->[0];
            my $startind=$#allpts;
            foreach my $thispt (@{$poly}[1..@{$poly}-1]) {
                push(@allsegs, $#allpts, $#allpts + 1);
                push(@allpts, $thispt);
                }
            if ($polycont->[0] eq 'polygon') { # add segment to close it
                push(@allsegs, $#allpts, $startind);                
                }

            }
        # add segments to C struct
        my $segs_added_count = $self->in->segmentlist(@allsegs);

        # Add region mark points and any attributes and area constraints to C struct
        if (@{$self->{poly}->{regions}}) {
            my $regions_added_count = $self->in->regionlist(map {grep defined, @{$_->[0]},$_->[1],$_->[2]} @{$self->{poly}->{regions}});
            }
        # Add hole mark points to C struct
        if (@{$self->{poly}->{holes}}) {
            my $holes_added_count = $self->in->holelist(map {@{$_}} @{$self->{poly}->{holes}});
            }
        }

    # add all points from PSLG, (possibly none)
    # and then any other points (not associated with segments)
    # into the C struct
    my $points_added_count = $self->in->pointlist(map {$_->[0],$_->[1]} (@allpts, @{$self->{poly}->{points}}));

    # set up attribute array if any points have more than 2 items (the coordinates) in list
    my $coords_plus_attrs = 2; # 2 for the coords - we'll skip over them when it's time
    foreach my $point (@allpts, @{$self->{poly}->{points}}) {
        if ($coords_plus_attrs < @{$point}) {$coords_plus_attrs = @{$point}}
        }
    if ($coords_plus_attrs > 2) {
        # Extend / fill in the attribute lists for any points 
        # that don't have the full set of attributes defined.
        # Set missing/undefined attributes to zero.
        foreach my $point (@allpts, @{$self->{poly}->{points}}) {
            if (@{$point} < $coords_plus_attrs) {
                foreach (2 .. $coords_plus_attrs - 1) {
                    if (!defined($point->[$_])) {$point->[$_]=0;}
                    }
                }
            }
        # put attributes into C struct
        $self->in->numberofpointattributes($coords_plus_attrs - 2);
        my $attributes_added_count = $self->in->pointattributelist(
            map {@{$_}[2 .. $coords_plus_attrs - 1]} (@allpts, @{$self->{poly}->{points}}));
        }

    # set up new triangulateio C structs to receive output
    $self->{out}    = new Math::Geometry::Delaunay::Triangulateio;
    $self->{vorout} = new Math::Geometry::Delaunay::Triangulateio;
    return;
    }


sub triangulate() {
    my $self = shift;
    my $dotopo = defined wantarray && !wantarray; # scalar or array return context 
    my $optstr = @_ ? join('',@_):'';
    $self->prepPoly($optstr);
    Math::Geometry::Delaunay::_triangulate($self->{optstr},$self->in->to_ptr,$self->out->to_ptr,$self->vorout->to_ptr);
    if (defined wantarray) { # else don't do expensive topology stuff if undef/void context
        if (wantarray && index($self->{optstr},'v') != -1) {
            return topology($self->out), topology($self->vorout);
            }
        return topology($self->out);
        }
    return;
    }

sub ltolol {($#_<$_[0])?():map [@_[$_*$_[0]+1..$_*$_[0]+$_[0]]],0..$#_/$_[0]-1}#perl.

sub nodes {
    my $self  = shift;
    my $fromVouty = @_ ? shift : 0;
    my $triio = $fromVouty ? $self->vorout : $self->out;
    my $cachetarg = $fromVouty ? 'voutnodes' : 'outnodes';
    if (@{$self->{poly}->{$cachetarg}} == 0) {
        my @nodeattributes;
        if ($triio->numberofpointattributes) {
            @nodeattributes = ltolol($triio->numberofpointattributes,$triio->pointattributelist);
            }
        @{$self->{poly}->{$cachetarg}} = ltolol(2,$triio->pointlist);
        for (my $i=0;$i<@nodeattributes;$i++) {
            push @{$self->{poly}->{$cachetarg}->[$i]}, @{$nodeattributes[$i]};
            } 
        if (!$fromVouty) {
            my @nodemarkers = $triio->pointmarkerlist;
            for (my $i=0;$i<@nodemarkers;$i++) {
                push @{$self->{poly}->{$cachetarg}->[$i]}, $nodemarkers[$i];
                }
            }
        }
    return $self->{poly}->{$cachetarg};
    }

sub elements {
    my $self = shift;
    my $triio = $self->out;
    my $nodes = $self->nodes;
    my @outelements;
    my @triangleattributes;
    if ($triio->numberoftriangleattributes) {
        @triangleattributes = ltolol($triio->numberoftriangleattributes,$triio->triangleattributelist);
        }
    @outelements = map {[map {$nodes->[$_]} @{$_}]} ltolol($triio->numberofcorners,$triio->trianglelist);
    for (my $i=0;$i<@triangleattributes;$i++) {
        push @{$outelements[$i]}, @{$triangleattributes[$i]};
        } 
    return \@outelements;
    }

sub segments {
    my $self = shift;
    my $triio = $self->out;
    my $nodes = $self->nodes;
    my @outsegments;
    my @segmentmarkers = $triio->segmentmarkerlist;
    @outsegments = map {[$nodes->[$_->[0]],$nodes->[$_->[1]]]} ltolol(2,$triio->segmentlist);
    for (my $i=0;$i<@segmentmarkers;$i++) {
        push @{$outsegments[$i]}, $segmentmarkers[$i];
        } 
    return \@outsegments;
    }
 
sub edges {
    my $self = shift;
    my $fromVouty = @_ ? shift : 0;
    my $triio = $fromVouty ? $self->vorout : $self->out;
    my $nodes = $self->nodes($fromVouty);
    my @outedges;
    @outedges = map {[map { $_==-1?0:$nodes->[$_]} @{$_}]} ltolol(2,$triio->edgelist);
    if (!$fromVouty) {
        my @edgemarkers = $triio->edgemarkerlist;
        for (my $i=0;$i<@edgemarkers;$i++) {
            push @{$outedges[$i]}, $edgemarkers[$i];
            }
        } 
    return \@outedges;
    }

sub vnodes {return $_[0]->nodes(1);}

sub vedges {
    my $self = shift;
    my $vedges = $self->edges(1);
    my $triio = $self->vorout;
    my @outrays;
    @outrays = ltolol(2,$triio->normlist);
    for (my $i=0;$i<@{$vedges};$i++) {
        # if one end was a ray (missing node ref)
        # look up the direction vector and use that as missing point
        # and set third element in edge array to true
        # as a flag to identify this edge as a ray 
        if (!$vedges->[$i]->[0]) {
            $vedges->[$i]->[0] = $outrays[$i];
            $vedges->[$i]->[2] = 1;
            }
        elsif (!$vedges->[$i]->[1]) {
            $vedges->[$i]->[1] = $outrays[$i];
            $vedges->[$i]->[2] = 2;
            }
        else {
            $vedges->[$i]->[2] = 0;
            }
        }
    return $vedges;
    }

sub topology {
    my $triio = shift;

    my $isVoronoi = 0; # we'll detect this when reading edges

    my @nodes = map {{ point => $_, attributes => [], marker => undef, elements => [], edges => [] , segments => []}} ltolol(2,$triio->pointlist);
    my @eles  = map {my $ele={ nodes=>[map {$nodes[$_]} @{$_}],edges => [], neighbors => [], marker => undef, attributes => [] };map {push @{$_->{elements}},$ele} @{$ele->{nodes}};$ele} ltolol($triio->numberofcorners,$triio->trianglelist);
    my $ecnt = 0; # The index for edges will be the link between Delaunay and Voronoi topology pairs.
    my @edges = map {my $edg={ nodes=>[map {$nodes[$_]} grep {$_>-1} @{$_}],marker => undef, elements => [], vector => undef, index => $ecnt++};foreach (@{$edg->{nodes}}) {push @{$_->{edges}   },$edg};$isVoronoi = 1 if ($_->[0] == -1 || $_->[1] == -1);$edg} ltolol(2,$triio->edgelist);
    my @segs  = map {my $edg={ nodes=>[map {$nodes[$_]}              @{$_}],marker => undef, elements => []                                   };foreach (@{$edg->{nodes}}) {push @{$_->{segments}},$edg};                                                   $edg} ltolol(2,$triio->segmentlist);
 
    my @elementattributes;
    if ($triio->numberoftriangleattributes) {
        @elementattributes = ltolol($triio->numberoftriangleattributes,$triio->triangleattributelist);
        }
    for (my $i=0;$i<@elementattributes;$i++) {
        $eles[$i]->{attributes} = $elementattributes[$i];
        }
    my @nodeattributes;
    if ($triio->numberofpointattributes) {
        @nodeattributes = ltolol($triio->numberofpointattributes,$triio->pointattributelist);
        }
    my @nodemarkers = $triio->pointmarkerlist; # always there for pslg, unlike attributes
    for (my $i=0;$i<@nodemarkers;$i++) {
        if ($triio->numberofpointattributes) {
            $nodes[$i]->{attributes} = $nodeattributes[$i];
            }
        $nodes[$i]->{marker} = $nodemarkers[$i];
        }
    my @edgemarkers = $triio->edgemarkerlist;
    for (my $i=0;$i<@edgemarkers;$i++) {
        $edges[$i]->{marker} = $edgemarkers[$i];
        }
    my @segmentmarkers = $triio->segmentmarkerlist; # because some can be internal to boundaries
    for (my $i=0;$i<@segmentmarkers;$i++) {
        $segs[$i]->{marker} = $segmentmarkers[$i];
        }
    my @neighs = ltolol(3,$triio->neighborlist);
    for (my $i=0;$i<@neighs;$i++) {
        $eles[$i]->{neighbors} = [map {$eles[$_]} grep {$_ != -1} @{$neighs[$i]}];
        }
    if ($isVoronoi) {
        my @edgevectors = ltolol(2,$triio->normlist); # voronoi ray vectors
        for (my $i=0;$i<@edgevectors;$i++) {
            $edges[$i]->{vector} = $edgevectors[$i];
            }
        }

    # cross reference elements and edges
    if (@edges) { # but only if edges were generated
        foreach my $ele (@eles) {
            for (my $i=-1;$i<@{$ele->{nodes}}-1;$i++) {
                foreach my $edge (@{$ele->{nodes}->[$i]->{edges}}) {
                    if ($ele->{nodes}->[$i+1] == $edge->{nodes}->[0] || $ele->{nodes}->[$i+1] == $edge->{nodes}->[1]) {
                        push @{$ele->{edges}}, $edge;
                        push @{$edge->{elements}}, $ele;
                        last;
                        }
                    }
                }
            }
        }
    return {
        nodes    => [@nodes],
        edges    => [@edges],
        segments => [@segs],
        elements => [@eles]
        };
    }

sub get_point_in_polygon {
    my $poly = shift;
    my $point_inside;

    my $bottom_left_index=0;
    my $maxy = $poly->[$bottom_left_index]->[1];
    for (my $i=1;$i<@{$poly};$i++) {
        if ($poly->[$i]->[1] <= $poly->[$bottom_left_index]->[1]) {
            if ($poly->[$i]->[1] < $poly->[$bottom_left_index]->[1] ||
                $poly->[$i]->[0] <  $poly->[$bottom_left_index]->[0]
               ) {
                $bottom_left_index = $i;
                }
            }
        if ($maxy < $poly->[$i]->[1]) { $maxy = $poly->[$i]->[1] }
        }
    my $prev_index = $bottom_left_index;
    my $next_index = -1 * @{$poly} + $bottom_left_index;
    --$prev_index 
        while ($poly->[$prev_index]->[0] == $poly->[$bottom_left_index]->[0] &&
               $poly->[$prev_index]->[1] == $poly->[$bottom_left_index]->[1]
              );
    ++$next_index
        while ($poly->[$next_index]->[0] == $poly->[$bottom_left_index]->[0] &&
               $poly->[$next_index]->[1] == $poly->[$bottom_left_index]->[1]
              );

    my @vec1 = ($poly->[$bottom_left_index]->[0] - $poly->[$prev_index]->[0],
                $poly->[$bottom_left_index]->[1] - $poly->[$prev_index]->[1]);
    my @vec2 = ($poly->[$next_index]->[0] - $poly->[$bottom_left_index]->[0],
                $poly->[$next_index]->[1] - $poly->[$bottom_left_index]->[1]);
    my $orient = (($vec1[0]*$vec2[1] - $vec2[0]*$vec1[1])>=0) ? 1:0;
    
    my $angle1 = atan2($poly->[$prev_index]->[1] - $poly->[$bottom_left_index]->[1], $poly->[$prev_index]->[0] - $poly->[$bottom_left_index]->[0]);
    my $angle2 = atan2($poly->[$next_index]->[1] - $poly->[$bottom_left_index]->[1], $poly->[$next_index]->[0] - $poly->[$bottom_left_index]->[0]);
    my $angle;        
    if ($orient) {$angle = $angle2 + ($angle1 - $angle2)/2;}
    else         {$angle = $angle1 + ($angle2 - $angle1)/2;}
    my $cosangle = cos($angle);
    my $sinangle = sin($angle);
    my $tanangle = $sinangle/$cosangle;

    my $adequate_distance = 1.1 * ($maxy - $poly->[$bottom_left_index]->[1]) *
                            ((abs($cosangle) < 0.00000001) ? 1 : 1/$sinangle);
    my $point_wayout = [$poly->[$bottom_left_index]->[0] + $adequate_distance * $cosangle, 
                        $poly->[$bottom_left_index]->[1] + $adequate_distance * $sinangle];
    my $ray = [$poly->[$bottom_left_index],$point_wayout];

    my @intersections = grep {
        # print "interdist: ",$_->[2],"\n";
        # Note the filtering out of very-near-zero distance intersections -
        # There should always be two intersections where the base of 
        # the test ray sits on the $poly->[$bottom_left_index] point:
        # one intersection for each adjacent segment. Better to get them
        # and ignore them, rather than try to move the base one way or the other
        # to avoid them.
        abs($_->[2]) > 10**-5} sort {$a->[2] <=> $b->[2]} seg_poly_intersections($ray,$poly,1);
    if (!@intersections) { 
        print "Warning: Failed to calculate hole or region indicator point."; 
        }
    elsif ($#intersections % 2 != 0) { 
        print "Warning: Calculated hole or region indicator point is not inside its polygon."; 
        }
    else {
        my $closest_intersection = $intersections[0];
        $point_inside = [($poly->[$bottom_left_index]->[0] + $closest_intersection->[0])/2,
                         ($poly->[$bottom_left_index]->[1] + $closest_intersection->[1])/2];
        #print "Calculated point inside: ",$point_inside->[0]," , ",$point_inside->[1],"\n";
        }
    return $point_inside, $ray;
    }

sub seg_poly_intersections {
    my $seg1 = shift;
    my $poly = shift;
    my $doDists = @_ ? shift : 0;
    my @intersections;
    for (my $i = -1; $i < @{$poly} - 1; $i++) {
        my $seg2 = [$poly->[$i],$poly->[$i+1]];

        my @segsegret;

        my $x1= $seg1->[0]->[0]; my $y1= $seg1->[0]->[1];
        my $x2= $seg1->[1]->[0]; my $y2= $seg1->[1]->[1];
        my $u1= $seg2->[0]->[0]; my $v1= $seg2->[0]->[1];
        my $u2= $seg2->[1]->[0]; my $v2= $seg2->[1]->[1];

        ##to maybe optimize for the case where segments are
        ##expected NOT to intersect most of the time
        #my @lowhix=($x2>$x1)?($x1,$x2):($x2,$x1);
        #my @lowhiu=($u2>$u1)?($u1,$u2):($u2,$u1);
        #if (
        #   $lowhix[0]>$lowhiu[1]
        #   ||
        #   $lowhix[1]<$lowhiu[0]   
        #   ) {
        #   return;
        #   }
        #my @lowhiy=($y2>$y1)?($y1,$y2):($y2,$y1);
        #my @lowhiv=($v2>$v1)?($v1,$v2):($v2,$v1);
        #if (
        #   $lowhiy[0]>$lowhiv[1]
        #   ||
        #   $lowhiy[1]<$lowhiv[0]   
        #   ) {
        #   return;
        #   }

        my $m1 = ($x2 eq $x1)?'Inf':($y2 - $y1)/($x2 - $x1);
        my $m2 = ($u2 eq $u1)?'Inf':($v2 - $v1)/($u2 - $u1);

        my $b1;
        my $b2;

        my  $xi;
        my $dm = $m1 - $m2;
        if    ($m1 eq 'Inf' && $m2 ne 'Inf') {$xi = $x1;$b2 = $v1 - ($m2 * $u1);}
        elsif ($m2 eq 'Inf' && $m1 ne 'Inf') {$xi = $u1;$b1 = $y1 - ($m1 * $x1);}
        elsif (abs($dm) > 0.000000000001) {
            $b1 = $y1 - ($m1 * $x1);
            $b2 = $v1 - ($m2 * $u1);    
            $xi=($b2-$b1)/$dm;
            }
        my @lowhiu=($u2>$u1)?($u1,$u2):($u2,$u1);
        if ($m1 ne 'Inf') {
            my @lowhix=($x2>$x1)?($x1,$x2):($x2,$x1);
            if ($m2 eq 'Inf' &&   ($u2<$lowhix[0] || $u2>$lowhix[1]) ) {
                #return ();
                }
            #if (!isNaN(xi) &&
            if (
                ($xi || $xi eq 0) &&
                ($xi < $lowhix[1] || $xi eq $lowhix[1]) && 
                ($xi > $lowhix[0] || $xi eq $lowhix[0]) &&
                ($xi < $lowhiu[1] || $xi eq $lowhiu[1]) && 
                ($xi > $lowhiu[0] || $xi eq $lowhiu[0])
                ) {
                my $y=($m1*$xi)+$b1;
                my @lowhiv=($v2>$v1)?($v1,$v2):($v2,$v1);
                if ($m2 eq 'Inf' &&
                    $y<$lowhiv[0] || $y>$lowhiv[1]
                    ) {
                    #return ();
                    }
                else {
                    push(@intersections,[$xi,$y]);
                    }
                }
            }
        elsif ($m2 ne 'Inf') {#so $m1 is Inf
            if ($x1 < $lowhiu[0] || $x1 > $lowhiu[1] && ! ($x1 eq $lowhiu[0] || $x1 eq $lowhiu[1])) {
                #return ();
                }
            my @lowhiy=($y2>$y1)?($y1,$y2):($y2,$y1);
            my @lowhiv=($v2>$v1)?($v1,$v2):($v2,$v1);
            my $yi = ($m2*$xi)+$b2;
            #print "$x1,$y1,$x2,$y2\n  $yi = ($m2*$xi)+$b2;\n";
            if (($yi || $yi eq 0) &&
                ($yi < $lowhiy[1] || $yi eq $lowhiy[1]) && 
                ($yi > $lowhiy[0] || $yi eq $lowhiy[0]) &&
                ($yi < $lowhiv[1] || $yi eq $lowhiv[1]) && 
                ($yi > $lowhiv[0] || $yi eq $lowhiv[0])
                ) {
                push(@intersections,[$xi,$yi]);
                }
            }
        }
    if ($doDists) {
        foreach my $int (@intersections) {
            push(@{$int}, sqrt( ($int->[0]-$seg1->[0]->[0])**2 + ($int->[1]-$seg1->[0]->[1])**2 ) );
            }
        }
    return @intersections;
    }


=head1 NAME

Math::Geometry::Delaunay - Quality Mesh Generator and Delaunay Triangulator

=head1 VERSION

Version 0.04

=cut

=head1 SYNOPSIS

=for html <div style="width:30%;float:right;display:inline-block;text-align:center;">
<svg viewBox="0 0 8 6" height="75px" preserveAspectRatio="xMinYMin meet" xmlns="http://www.w3.org/2000/svg" version="1.1">
<g transform="scale(1,-1) translate(0,-6)">
<circle r="0.15" cx="1" cy="1" style="fill:blue"/>
<circle r="0.15" cx="7" cy="1" style="fill:blue"/>
<circle r="0.15" cx="7" cy="3" style="fill:blue"/>
<circle r="0.15" cx="3" cy="3" style="fill:blue"/>
<circle r="0.15" cx="3" cy="5" style="fill:blue"/>
<circle r="0.15" cx="1" cy="5" style="fill:blue"/>
</g></svg><br/>
<br/><small>input vertices</small><br/>
<svg viewBox="0 0 8 6" height="75px" preservAspectRatio="minXminY meet" xmlns="http://www.w3.org/2000/svg" version="1.1">
<g transform="scale(1,-1) translate(0,-6)">
<line x1="1" y1="1" x2="3" y2="3" style="stroke:gray;stroke-width:0.1;"/>
<line x1="3" y1="3" x2="1" y2="5" style="stroke:gray;stroke-width:0.1;"/>
<line x1="1" y1="5" x2="1" y2="1" style="stroke:gray;stroke-width:0.1;"/>
<line x1="1" y1="1" x2="7" y2="1" style="stroke:gray;stroke-width:0.1;"/>
<line x1="7" y1="1" x2="3" y2="3" style="stroke:gray;stroke-width:0.1;"/>
<line x1="3" y1="3" x2="3" y2="5" style="stroke:gray;stroke-width:0.1;"/>
<line x1="3" y1="5" x2="1" y2="5" style="stroke:gray;stroke-width:0.1;"/>
<line x1="3" y1="3" x2="7" y2="3" style="stroke:gray;stroke-width:0.1;"/>
<line x1="7" y1="3" x2="3" y2="5" style="stroke:gray;stroke-width:0.1;"/>
<line x1="7" y1="1" x2="7" y2="3" style="stroke:gray;stroke-width:0.1;"/>
<circle cx="1" cy="1" r="0.150" style="fill:blue;"/>
<circle cx="7" cy="1" r="0.150" style="fill:blue;"/>
<circle cx="7" cy="3" r="0.150" style="fill:blue;"/>
<circle cx="3" cy="3" r="0.150" style="fill:blue;"/>
<circle cx="3" cy="5" r="0.150" style="fill:blue;"/>
<circle cx="1" cy="5" r="0.150" style="fill:blue;"/>
</g></svg>
<br/><small>Delaunay triangulation</small><br/>
<svg viewBox="0 0 8 6" height="75px" preservAspectRatio="minXminY meet" xmlns="http://www.w3.org/2000/svg" xmlns:svg="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:slic3r="http://slic3r.org/namespaces/slic3r">
<g transform="scale(1,-1) translate(0,-5)">
<line x1="1" y1="3" x2="4" y2="0" style="stroke:gray;stroke-width:0.1;"/>
<line x1="1" y1="3" x2="2" y2="4" style="stroke:gray;stroke-width:0.1;"/>
<line x1="4" y1="0" x2="5" y2="2" style="stroke:gray;stroke-width:0.1;"/>
<line x1="2" y1="4" x2="5" y2="4" style="stroke:gray;stroke-width:0.1;"/>
<line x1="5" y1="4" x2="5" y2="2" style="stroke:gray;stroke-width:0.1;"/>
<line x1="1" y1="3" x2="-139" y2="3" style="stroke:gray;stroke-width:0.1;"/>
<line x1="4" y1="0" x2="4" y2="-210" style="stroke:gray;stroke-width:0.1;"/>
<line x1="2" y1="4" x2="2" y2="74" style="stroke:gray;stroke-width:0.1;"/>
<line x1="5" y1="4" x2="75" y2="144" style="stroke:gray;stroke-width:0.1;"/>
<line x1="5" y1="2" x2="75" y2="2" style="stroke:gray;stroke-width:0.1;"/>
<circle cx="1" cy="3" r="0.15" style="fill:blue;"/>
<circle cx="4" cy="0" r="0.15" style="fill:blue;"/>
<circle cx="2" cy="4" r="0.15" style="fill:blue;"/>
<circle cx="5" cy="4" r="0.15" style="fill:blue;"/>
<circle cx="5" cy="2" r="0.15" style="fill:blue;"/>
</g></svg>
<br/><small>Voronoi diagram</small>
</div>

    use Math::Geometry::Delaunay qw(TRI_CCDT);

    # generate Delaunay triangulation
    # and Voronoi diagram for a point set

    my $point_set = [ [1,1], [7,1], [7,3],
                      [3,3], [3,5], [1,5] ];

    my $tri = new Math::Geometery::Delaunay('ev');
    $tri->addPoints($point_set);
    
    # called in void context
    $tri->triangulate();
    # populates the following lists

    $tri->elements(); # triangles
    $tri->nodes();    # points
    $tri->edges();    # triangle edges
    $tri->vnodes();   # Voronoi diagram points
    $tri->vedges();   # Voronoi edges and rays


=for html <br clear="all"/><div style="margin-top:30px;width:30%;float:right;display:inline-block;text-align:center;">
<svg viewBox="0 0 8 6" height="75px" preserveAspectRatio="xMinYMin meet" xmlns="http://www.w3.org/2000/svg" version="1.1">
<g transform="scale(1,-1) translate(0,-6)">
<line x1="7" y1="1" x2="1" y2="1" style="stroke:blue;stroke-width:0.1;"/>
<line x1="7" y1="3" x2="7" y2="1" style="stroke:blue;stroke-width:0.1;"/>
<line x1="3" y1="3" x2="7" y2="3" style="stroke:blue;stroke-width:0.1;"/>
<line x1="3" y1="3" x2="3" y2="5" style="stroke:blue;stroke-width:0.1;"/>
<line x1="1" y1="5" x2="3" y2="5" style="stroke:blue;stroke-width:0.1;"/>
<line x1="1" y1="1" x2="1" y2="5" style="stroke:blue;stroke-width:0.1;"/>
<circle cx="1" cy="1" r="0.15" style="fill:black;"/>
<circle cx="7" cy="1" r="0.15" style="fill:black;"/>
<circle cx="7" cy="3" r="0.15" style="fill:black;"/>
<circle cx="3" cy="3" r="0.15" style="fill:black;"/>
<circle cx="3" cy="5" r="0.15" style="fill:black;"/>
<circle cx="1" cy="5" r="0.15" style="fill:black;"/>
</g></svg>
<br/><small>input PSLG</small><br/>
<svg viewBox="0 0 8 6" height="75px" preserveAspectRatio="xMinYMin meet" style="margin-top:23px;" xmlns="http://www.w3.org/2000/svg" version="1.1">
<g transform="scale(1,-1) translate(0,-6)">
<line x1="1.55" y1="2.5" x2="1" y2="3" style="stroke:gray;stroke-width:0.1;"/>
<line x1="1" y1="3" x2="1" y2="2" style="stroke:gray;stroke-width:0.1;"/>
<line x1="1" y1="2" x2="1.55" y2="2.5" style="stroke:gray;stroke-width:0.1;"/>
<line x1="4.75" y1="1.72322767215881" x2="4" y2="1" style="stroke:gray;stroke-width:0.1;"/>
<line x1="4" y1="1" x2="5.5" y2="1" style="stroke:gray;stroke-width:0.1;"/>
<line x1="5.5" y1="1" x2="4.75" y2="1.72322767215881" style="stroke:gray;stroke-width:0.1;"/>
<line x1="3" y1="5" x2="2" y2="5" style="stroke:gray;stroke-width:0.1;"/>
<line x1="2" y1="5" x2="2.25" y2="4.25" style="stroke:gray;stroke-width:0.1;"/>
<line x1="2.25" y1="4.25" x2="3" y2="5" style="stroke:gray;stroke-width:0.1;"/>
<line x1="3" y1="4" x2="2.25" y2="4.25" style="stroke:gray;stroke-width:0.1;"/>
<line x1="2.25" y1="4.25" x2="2" y2="3.5" style="stroke:gray;stroke-width:0.1;"/>
<line x1="2" y1="3.5" x2="3" y2="4" style="stroke:gray;stroke-width:0.1;"/>
<line x1="2.25" y1="2.75" x2="3" y2="3" style="stroke:gray;stroke-width:0.1;"/>
<line x1="3" y1="3" x2="2" y2="3.5" style="stroke:gray;stroke-width:0.1;"/>
<line x1="2" y1="3.5" x2="2.25" y2="2.75" style="stroke:gray;stroke-width:0.1;"/>
<line x1="2.93181818181818" y1="1.95454545454545" x2="3.71892953370153" y2="1.8730674509337" style="stroke:gray;stroke-width:0.1;"/>
<line x1="3.71892953370153" y1="1.8730674509337" x2="3.5" y2="2.52618854085072" style="stroke:gray;stroke-width:0.1;"/>
<line x1="3.5" y1="2.52618854085072" x2="2.93181818181818" y2="1.95454545454545" style="stroke:gray;stroke-width:0.1;"/>
<line x1="3" y1="3" x2="3" y2="4" style="stroke:gray;stroke-width:0.1;"/>
<line x1="2.25" y1="4.25" x2="1.33333333333333" y2="4.13888888888889" style="stroke:gray;stroke-width:0.1;"/>
<line x1="1.33333333333333" y1="4.13888888888889" x2="2" y2="3.5" style="stroke:gray;stroke-width:0.1;"/>
<line x1="3" y1="4" x2="3" y2="5" style="stroke:gray;stroke-width:0.1;"/>
<line x1="2" y1="5" x2="1.33333333333333" y2="4.13888888888889" style="stroke:gray;stroke-width:0.1;"/>
<line x1="2" y1="5" x2="1" y2="5" style="stroke:gray;stroke-width:0.1;"/>
<line x1="1" y1="5" x2="1.33333333333333" y2="4.13888888888889" style="stroke:gray;stroke-width:0.1;"/>
<line x1="1" y1="5" x2="1" y2="3" style="stroke:gray;stroke-width:0.1;"/>
<line x1="1" y1="3" x2="1.33333333333333" y2="4.13888888888889" style="stroke:gray;stroke-width:0.1;"/>
<line x1="1" y1="3" x2="2" y2="3.5" style="stroke:gray;stroke-width:0.1;"/>
<line x1="2.93181818181818" y1="1.95454545454545" x2="1.8375" y2="1.63125" style="stroke:gray;stroke-width:0.1;"/>
<line x1="1.8375" y1="1.63125" x2="2.5" y2="1" style="stroke:gray;stroke-width:0.1;"/>
<line x1="2.5" y1="1" x2="2.93181818181818" y2="1.95454545454545" style="stroke:gray;stroke-width:0.1;"/>
<line x1="2.25" y1="2.75" x2="1.8375" y2="1.63125" style="stroke:gray;stroke-width:0.1;"/>
<line x1="2.93181818181818" y1="1.95454545454545" x2="2.25" y2="2.75" style="stroke:gray;stroke-width:0.1;"/>
<line x1="2.5" y1="1" x2="3.25" y2="1.23566017316017" style="stroke:gray;stroke-width:0.1;"/>
<line x1="3.25" y1="1.23566017316017" x2="2.93181818181818" y2="1.95454545454545" style="stroke:gray;stroke-width:0.1;"/>
<line x1="1.8375" y1="1.63125" x2="1" y2="1" style="stroke:gray;stroke-width:0.1;"/>
<line x1="1" y1="1" x2="2.5" y2="1" style="stroke:gray;stroke-width:0.1;"/>
<line x1="1.55" y1="2.5" x2="2" y2="3.5" style="stroke:gray;stroke-width:0.1;"/>
<line x1="2.25" y1="2.75" x2="1.55" y2="2.5" style="stroke:gray;stroke-width:0.1;"/>
<line x1="1.55" y1="2.5" x2="1.8375" y2="1.63125" style="stroke:gray;stroke-width:0.1;"/>
<line x1="1.8375" y1="1.63125" x2="1" y2="2" style="stroke:gray;stroke-width:0.1;"/>
<line x1="1" y1="2" x2="1" y2="1" style="stroke:gray;stroke-width:0.1;"/>
<line x1="2.5" y1="1" x2="4" y2="1" style="stroke:gray;stroke-width:0.1;"/>
<line x1="4" y1="1" x2="3.25" y2="1.23566017316017" style="stroke:gray;stroke-width:0.1;"/>
<line x1="2.93181818181818" y1="1.95454545454545" x2="3" y2="3" style="stroke:gray;stroke-width:0.1;"/>
<line x1="3.71892953370153" y1="1.8730674509337" x2="4.75" y2="1.72322767215881" style="stroke:gray;stroke-width:0.1;"/>
<line x1="4.75" y1="1.72322767215881" x2="4.5" y2="2.43504117935408" style="stroke:gray;stroke-width:0.1;"/>
<line x1="4.5" y1="2.43504117935408" x2="3.71892953370153" y2="1.8730674509337" style="stroke:gray;stroke-width:0.1;"/>
<line x1="7" y1="2" x2="6.5096079838992" y2="2.5" style="stroke:gray;stroke-width:0.1;"/>
<line x1="6.5096079838992" y1="2.5" x2="6.36316415068211" y2="1.5" style="stroke:gray;stroke-width:0.1;"/>
<line x1="6.36316415068211" y1="1.5" x2="7" y2="2" style="stroke:gray;stroke-width:0.1;"/>
<line x1="4" y1="3" x2="3" y2="3" style="stroke:gray;stroke-width:0.1;"/>
<line x1="3" y1="3" x2="3.5" y2="2.52618854085072" style="stroke:gray;stroke-width:0.1;"/>
<line x1="3.5" y1="2.52618854085072" x2="4" y2="3" style="stroke:gray;stroke-width:0.1;"/>
<line x1="3.25" y1="1.23566017316017" x2="3.71892953370153" y2="1.8730674509337" style="stroke:gray;stroke-width:0.1;"/>
<line x1="3.71892953370153" y1="1.8730674509337" x2="4" y2="1" style="stroke:gray;stroke-width:0.1;"/>
<line x1="5" y1="3" x2="4.5" y2="2.43504117935408" style="stroke:gray;stroke-width:0.1;"/>
<line x1="4.5" y1="2.43504117935408" x2="5.42596576190246" y2="1.67372070106296" style="stroke:gray;stroke-width:0.1;"/>
<line x1="5.42596576190246" y1="1.67372070106296" x2="5" y2="3" style="stroke:gray;stroke-width:0.1;"/>
<line x1="5.42596576190246" y1="1.67372070106296" x2="4.75" y2="1.72322767215881" style="stroke:gray;stroke-width:0.1;"/>
<line x1="5.5" y1="1" x2="5.42596576190246" y2="1.67372070106296" style="stroke:gray;stroke-width:0.1;"/>
<line x1="4" y1="3" x2="4.5" y2="2.43504117935408" style="stroke:gray;stroke-width:0.1;"/>
<line x1="5" y1="3" x2="4" y2="3" style="stroke:gray;stroke-width:0.1;"/>
<line x1="4.5" y1="2.43504117935408" x2="3.5" y2="2.52618854085072" style="stroke:gray;stroke-width:0.1;"/>
<line x1="5.5" y1="1" x2="6.36316415068211" y2="1.5" style="stroke:gray;stroke-width:0.1;"/>
<line x1="6.36316415068211" y1="1.5" x2="5.42596576190246" y2="1.67372070106296" style="stroke:gray;stroke-width:0.1;"/>
<line x1="5" y1="3" x2="5.89643888151886" y2="2.16160972037971" style="stroke:gray;stroke-width:0.1;"/>
<line x1="5.89643888151886" y1="2.16160972037971" x2="6.5096079838992" y2="2.5" style="stroke:gray;stroke-width:0.1;"/>
<line x1="6.5096079838992" y1="2.5" x2="5" y2="3" style="stroke:gray;stroke-width:0.1;"/>
<line x1="7" y1="1" x2="7" y2="2" style="stroke:gray;stroke-width:0.1;"/>
<line x1="6.36316415068211" y1="1.5" x2="7" y2="1" style="stroke:gray;stroke-width:0.1;"/>
<line x1="5.42596576190246" y1="1.67372070106296" x2="5.89643888151886" y2="2.16160972037971" style="stroke:gray;stroke-width:0.1;"/>
<line x1="6.5096079838992" y1="2.5" x2="7" y2="3" style="stroke:gray;stroke-width:0.1;"/>
<line x1="7" y1="3" x2="5" y2="3" style="stroke:gray;stroke-width:0.1;"/>
<line x1="7" y1="2" x2="7" y2="3" style="stroke:gray;stroke-width:0.1;"/>
<line x1="5.89643888151886" y1="2.16160972037971" x2="6.36316415068211" y2="1.5" style="stroke:gray;stroke-width:0.1;"/>
<line x1="5.5" y1="1" x2="7" y2="1" style="stroke:gray;stroke-width:0.1;"/>
<line x1="7" y1="1" x2="5.5" y2="1" style="stroke:blue;stroke-width:0.1;"/>
<line x1="7" y1="3" x2="7" y2="2" style="stroke:blue;stroke-width:0.1;"/>
<line x1="3" y1="3" x2="4" y2="3" style="stroke:blue;stroke-width:0.1;"/>
<line x1="3" y1="4" x2="3" y2="5" style="stroke:blue;stroke-width:0.1;"/>
<line x1="1" y1="5" x2="2" y2="5" style="stroke:blue;stroke-width:0.1;"/>
<line x1="1" y1="1" x2="1" y2="2" style="stroke:blue;stroke-width:0.1;"/>
<line x1="3" y1="4" x2="3" y2="3" style="stroke:blue;stroke-width:0.1;"/>
<line x1="2" y1="5" x2="3" y2="5" style="stroke:blue;stroke-width:0.1;"/>
<line x1="1" y1="3" x2="1" y2="5" style="stroke:blue;stroke-width:0.1;"/>
<line x1="1" y1="2" x2="1" y2="3" style="stroke:blue;stroke-width:0.1;"/>
<line x1="4" y1="1" x2="2.5" y2="1" style="stroke:blue;stroke-width:0.1;"/>
<line x1="2.5" y1="1" x2="1" y2="1" style="stroke:blue;stroke-width:0.1;"/>
<line x1="5" y1="3" x2="7" y2="3" style="stroke:blue;stroke-width:0.1;"/>
<line x1="4" y1="3" x2="5" y2="3" style="stroke:blue;stroke-width:0.1;"/>
<line x1="5.5" y1="1" x2="4" y2="1" style="stroke:blue;stroke-width:0.1;"/>
<line x1="7" y1="2" x2="7" y2="1" style="stroke:blue;stroke-width:0.1;"/>
<circle cx="1" cy="1" r="0.15" style="fill:black;"/>
<circle cx="7" cy="1" r="0.15" style="fill:black;"/>
<circle cx="7" cy="3" r="0.15" style="fill:black;"/>
<circle cx="3" cy="3" r="0.15" style="fill:black;"/>
<circle cx="3" cy="5" r="0.15" style="fill:black;"/>
<circle cx="1" cy="5" r="0.15" style="fill:black;"/>
<circle cx="1" cy="3" r="0.15" style="fill:black;"/>
<circle cx="3" cy="4" r="0.15" style="fill:black;"/>
<circle cx="2" cy="5" r="0.15" style="fill:black;"/>
<circle cx="2" cy="3.5" r="0.15" style="fill:black;"/>
<circle cx="2.25" cy="4.25" r="0.15" style="fill:black;"/>
<circle cx="1.33" cy="4.138" r="0.15" style="fill:black;"/>
<circle cx="2.9318" cy="1.95" r="0.15" style="fill:black;"/>
<circle cx="4" cy="1" r="0.15" style="fill:black;"/>
<circle cx="2.25" cy="2.75" r="0.15" style="fill:black;"/>
<circle cx="1" cy="2" r="0.15" style="fill:black;"/>
<circle cx="1.55" cy="2.5" r="0.15" style="fill:black;"/>
<circle cx="1.8375" cy="1.63125" r="0.15" style="fill:black;"/>
<circle cx="2.5" cy="1" r="0.15" style="fill:black;"/>
<circle cx="3.25" cy="1.235" r="0.15" style="fill:black;"/>
<circle cx="5.5" cy="1" r="0.15" style="fill:black;"/>
<circle cx="4" cy="3" r="0.15" style="fill:black;"/>
<circle cx="5" cy="3" r="0.15" style="fill:black;"/>
<circle cx="3.7189" cy="1.873" r="0.15" style="fill:black;"/>
<circle cx="4.75" cy="1.723" r="0.15" style="fill:black;"/>
<circle cx="3.5" cy="2.526" r="0.15" style="fill:black;"/>
<circle cx="5.896" cy="2.1616" r="0.15" style="fill:black;"/>
<circle cx="4.5" cy="2.435" r="0.15" style="fill:black;"/>
<circle cx="5.4259" cy="1.6737" r="0.15" style="fill:black;"/>
<circle cx="7" cy="2" r="0.15" style="fill:black;"/>
<circle cx="6.5096" cy="2.5" r="0.15" style="fill:black;"/>
<circle cx="6.36316" cy="1.5" r="0.15" style="fill:black;"/>
</g></svg>
<br/><small>output mesh</small><br/>
<svg viewBox="0 0 8 6" height="75px" preserveAspectRatio="xMinYMin meet" style="margin-top:23px;" xmlns="http://www.w3.org/2000/svg" version="1.1">
<g transform="scale(1,-1) translate(0,-6)">
<g style="fill:black">
<polygon points="1.55,2.5,1,3,1,2"/>
<polygon points="4.75,1.723,4,1,5.5,1"/>
<polygon points="3,5,2,5,2.25,4.25"/>
<polygon points="3,3,3,4,2,3.5"/>
<polygon points="3,5,2.25,4.25,3,4"/>
<polygon points="1.333,4.139,2,5,1,5"/>
<polygon points="1,5,1,3,1.333,4.139"/>
<polygon points="1.837,1.631,1,1,2.5,1"/>
<polygon points="1,1,1.837,1.631,1,2"/>
<polygon points="2.5,1,4,1,3.25,1.236"/>
<polygon points="4,3,3,3,3.5,2.526"/>
<polygon points="4,3,4.5,2.435,5,3"/>
<polygon points="7,1,7,2,6.363,1.5"/>
<polygon points="5,3,6.51,2.5,7,3"/>
<polygon points="7,3,6.51,2.5,7,2"/>
<polygon points="5.5,1,7,1,6.363,1.5"/>
</g>
<g style="fill:gray">
<polygon points="3,4,2.25,4.25,2,3.5"/>
<polygon points="2.25,2.75,3,3,2,3.5"/>
<polygon points="2,5,1.333,4.139,2.25,4.25"/>
<polygon points="1,3,2,3.5,1.333,4.139"/>
<polygon points="2.932,1.955,1.837,1.631,2.5,1"/>
<polygon points="2.5,1,3.25,1.236,2.932,1.955"/>
<polygon points="1.55,2.5,2,3.5,1,3"/>
<polygon points="2.25,2.75,2.932,1.955,3,3"/>
<polygon points="1.55,2.5,1,2,1.837,1.631"/>
<polygon points="7,2,6.51,2.5,6.363,1.5"/>
<polygon points="4,1,4.75,1.723,3.719,1.873"/>
<polygon points="4,1,3.719,1.873,3.25,1.236"/>
<polygon points="5,3,4.5,2.435,5.426,1.674"/>
<polygon points="5.426,1.674,4.75,1.723,5.5,1"/>
<polygon points="2.932,1.955,3.5,2.526,3,3"/>
<polygon points="5.5,1,6.363,1.5,5.426,1.674"/>
<polygon points="5,3,5.896,2.162,6.51,2.5"/>
<polygon points="4.5,2.435,4,3,3.5,2.526"/>
<polygon points="5.426,1.674,5.896,2.162,5,3"/>
</g>
</g></svg>
<br/><small>something interesting<br/>extracted from topology</small>
</div>

    # quality mesh of a planar straight line graph
    # with cross-referenced topological output

    my $tri = new Math::Geometery::Delaunay('e');
    $tri->addPolygon($point_set);
    $tri->minimum_angle(23);

    # called in scalar context
    my $mesh_topology = $tri->triangulate(TRI_CCDT);
    # returns cross-referenced topology

    # make two lists of triangles that touch boundary segments

    my @tris_with_boundary_segment;
    my @tris_with_boundary_point;

    foreach my $triangle (@{$mesh_topology->{elements}}) {
            my $nodes_on_boundary_count = ( 
                grep $_->{marker} == 1,
                @{$triangle->{nodes}} 
                );
            if ($nodes_on_boundary_count == 2) {
                push @tris_with_boundary_segment, $triangle;
                }
            elsif ($nodes_on_boundary_count == 1) {
                push @tris_with_boundary_point, $triangle;
                }
            }
            

=for html <br clear="all"/>

=head1 DESCRIPTION

This is a Perl interface to the Jonathan Shewchuk's Triangle library.

"Triangle generates exact Delaunay triangulations, constrained Delaunay 
triangulations, conforming Delaunay triangulations, Voronoi diagrams, and 
high-quality triangular meshes. The latter can be generated with no small or 
large angles, and are thus suitable for finite element analysis." 
-- from L<http://www.cs.cmu.edu/~quake/triangle.html>

=head1 EXPORTS

Triangle has several option switches that can be used in different combinations 
to choose a class of triangulation and then configure options within that class.
To clarify the composition of option strings, or just to give you a head start, 
a few constants are supplied to configure different classes of mesh output.

    TRI_CONSTRAINED  = 'Y'    for "Constrained Delaunay"
    TRI_CONFORMING   = 'Dq0'  for "Conforming Delaunay"
    TRI_CCDT         = 'q'    for "Constrained Conforming Delaunay"
    TRI_VORONOI      = 'v'    to generate the Voronoi diagram

For an illustration of these terms, see: 
L<http://www.cs.cmu.edu/~quake/triangle.defs.html>

=head1 CONSTRUCTOR

=head2 new

The constructor returns a Math::Geometry::Delaunay object, and optionally takes 
a string, or list of strings, comprised of switches corresponding to Triangle's 
command line switches, documented here: 
L<http://www.cs.cmu.edu/~quake/triangle.switch.html>

    my $tri = Math::Geometry::Delaunay->new();

    my $tri = Math::Geometry::Delaunay->new('pzq0eQ');

    my $tri = Math::Geometry::Delaunay->new(TRI_CCDT, 'q15', 'a3.5');


Options set by switches passed to C<new()> may be overridden later by the 
corresponding option-setting methods, any time before C<triangulate()> 
is invoked.

=head1 MESH GENERATION

=head2 triangulate

Run the triangulation and either populate the object's output lists, or return 
a hash reference giving access to a cross-referenced representation of the mesh 
topology.

=head3 list output

After triangulate is invoked in void context, the output mesh data can be 
retrieved from the following methods, all of which return a reference to an 
array.

    $tri->triangulate(); # void context - no return value requested
    # output lists now available
    $points  = $tri->nodes();    # array of vertices
    $tris    = $tri->elements(); # array of triangles
    $edges   = $tri->edges();    # all the triangle edges
    $segs    = $tri->segments(); # the PSLG segments
    $vpoints = $tri->vnodes();   # points in the voronoi diagram
    $vedges  = $tri->vedges();   # edges in the voronoi diagram

Data may not be available for all lists, depending on which option switches were
used. By default, nodes and elements are generated, while edges are not.

The members of the lists returned have these formats:

    nodes:    [x, y, < zero or more attributes >, < boundary marker >]

    elements: [[x0, y0], [x1, y1], [x2, y2],
                < another three vertices, if "o2" switch used >,
                < zero or more attributes >
                ]
    edges:    [[x0, y0], [x1, y1], < boundary marker >]

    segments: [[x0, y0], [x1, y1], < boundary marker >]

    vnodes:   [x, y, < zero or more attributes >]

    vedges:   [< vertex or vector >, < vertex or vector >, < ray flag >]

Boundary markers are 1 or 0. An edge or segment with only one end on a boundary 
has boundary marker 0.

The ray flag is 0 if the edge is not a ray, or 1 or 2, to indicate 
which vertex is actually a unit vector indicating the direction of the ray.

Import of the mesh data from the C data structures will be deferred until
actually requested from the list fetching methods above. For speed and 
lower memory footprint, access only what you need, and consider suppressing 
output you don't need with option switches.

=head3 topological output

When triangulate is invoked in scalar or array context, it returns a hash ref 
containing the cross-referenced nodes, elements, edges, and PSLG segments of the
triangulation. In array context, with the "v" switch enabled, the Voronoi
topology is the second item returned.

    my $topology = $tri->triangulate();

    $topology now looks like this:
    
    {
    nodes    => [
                  { # a node
                  point      => [x0, x1],
                  edges      => [edgeref, ...],
                  segments   => [edgeref, ...], # a subset of edges
                  elements   => [elementref, ...],
                  marker     => 1 or 0 or undefined, # boundary marker
                  attributes => [attr0, ...]
                  },
                  ... more nodes like that
                    
                ],
    elements => [
                  { # a triangle
                  nodes      => [noderef0, noderef1, noderef2],
                  edges      => [edgeref0, edgeref1],
                  neighbors  => [neighref0, neighref1, neighref2],
                  attributes => [attrib0, ...]
                  },
                  ... more triangles like that
                ],
    edges    => [
                  {
                  nodes    => [noderef0, noderef1], # only one for a ray
                  elements => [elemref0, elemref1], # one if on boundary
                  vector   => undefined or [x, y],  # ray direction
                  marker   => 1 or 0 or undefined,  # boundary marker
                  index    => <integer> # edge's index in edge list
                  },
                  ... more edges like that
                    
                ],
    segments => [
                  {
                  nodes    => [noderef0, noderef1],
                  elements => [elemref0, elemref1], # one if on boundary
                  marker   => 1 or 0 or undefined   # boundary marker
                  },
                  ... more segments
                ]
    }

=head3 cross-referencing Delaunay and Voronoi

Corresponding edges in the Delaunay and Voronoi outputs have the same index
number in their respective edge lists. 

In the topological output, any edge in a triangulation has a record of its own 
index number that can by used to look up the corresponding edge in the Voronoi 
diagram topology, or vice versa, like so:

    ($topo, $voronoi_topo) = $tri->triangulate('ev');
    
    # get an edge reference where it's not obvious what the edge's index is
    
    $delaunay_edge = $topo->{nodes}->[-1]->{edges}->[-1];
    
    # this gets a reference to the corresponding edge in the Voronoi diagram
    
    $voronoi_edge = $voronoi_topo->{edges}->[$delaunay_edge->{index}];

=head1 METHODS TO SET SOME Triangle OPTIONS

=head2 area_constraint

Corresponds to the "a" switch.

With one argument, sets the maximum triangle area constraint for the 
triangulation. Returns the value supplied. With no argument, returns the 
current area constraint.

Passing -1 to C<area_constraint()> will disable the global area constraint.

=head2 minimum_angle

Corresponds to the "q" switch.

With one argument, sets the minimum angle allowed for triangles added in the
triangulation. Returns the value supplied. With no argument, returns the
current minimum angle constraint.

Passing -1 to C<minimum_angle()> will cause the "q" switch to be omitted from
the option string.

=head2 doEdges, doVoronoi, doNeighbors

These methods simply add or remove the corresponding letters from the 
option string. Pass in a true or false value to enable or disable.
Invoke with no argument to read the current state.

=head2 quiet, verbose

Triangle prints a basic summary of the meshing operation to STDOUT unless
the "Q" switch is present. This module includes the "Q" switch by default, but
you can override this by passing a false value to C<quiet()>.

If you would like to see even more output regarding the triangulation process,
there are are three levels of verbosity configurable with repeated "V"
switches. Passing a number from 1 to 3 to the C<verbose()> method will enable 
the corresponding level of verbosity.

=head1 METHODS TO ADD VERTICES AND SEGMENTS

=head2 addVertices, addPoints

Takes a reference to an array of vertices, each vertex itself an reference to
an array containing two coordinates and zero or more attributes. Attributes
are floating point numbers.
    
    # vertex format
    # [x, y, < zero or more attributes as floating point numbers >]

    $tri->addPoints([[$x0, $y0], [$x1, $y1], ... ]);

Use addVertices to add vertices that are not part of a PSLG. 
Use addPoints to add points that are not part of a polygon or polyline.
In other words, they do the same thing.

=head2 addSegments

Takes a reference to an array of segments.

    # segment format
    # [[$x0, $y0], [$x1, $y1]]

    $tri->addSegments([ $segment0, $segment1, ... ]);

If your segments are contiguous, it's better to use addPolyline, or addPolygon.

This method is provided because some point and polygon processing algorithms
result in segments that represent polygons, but list the segments in a 
non-contiguous order, and have shared vertices repeated in each segment's record.

The segments added with this method will be checked for duplicate vertices, and 
references to these will be merged.

Triangle can handle duplicate vertices, but we would rather not feed them in on 
purpose.

=head2 addPolyline

Takes a reference to an array of vertices describing a curve. 
Creates PSLG segments for each pair of adjacent vertices. Adds the
new segments and vertices to the triangulation input.

    $tri->addPolyline([$vertex0, $vertex1, $vertex2, ...]);

=head2 addPolygon

Takes a reference to an array of vertices describing a polygon. 
Creates PSLG segments for each pair of adjacent vertices
and creates and additional segment linking the last vertex to
the first,to close the polygon.  Adds the new segments and vertices 
to the triangulation input.

    $tri->addPolygon([$vertex0, $vertex1, $vertex2, ...]);

=head2 addHole

Like addPolygon, but describing a hole or concavity - an area of the output mesh
that should not be triangulated. 

There are two ways to specify a hole. Either provide a list of vertices, like
for addPolygon, or provide a single vertex that lies inside of a polygon, to
identify that polygon as a hole.

    # first way
    $tri->addHole([$vertex0, $vertex1, $vertex2, ...]);

    # second way
    $tri->addPolygon( [ [0,0], [1,0], [1,1], [0,1] ] );
    $tri->addHole( [0.5,0.5] );

Hole marker points can also be used, in combination with the "c" option, to
cause or preserve concavities in a boundary when Triangle would otherwise
enclose a PSLG in a convex hull.

=head2 addRegion

Takes a polygon describing a region, and an attribute or area constraint. With
both the "A" and "a" switches in effect, three arguments allow you to specify
both an attribute and an optional area constraint.

The first argument may alternately be a single vertex that lies inside of 
another polygon, to identify that polygon as a region.

To be used in conjunction with the "A" and "a" switches.

    # with the "A" switch
    $tri->addRegion(\@polygon, < attribute > );
    
    # with the "a" switch
    $tri->addRegion(\@polygon, < area constraint > );

    # with both "Aa"
    $tri->addRegion(\@polygon, < attribute >, < area constraint > );

If the "A" switch is used, each triangle generated within the bounds of a region
will have that region's attribute added to the end of the triangle's 
attributes list, while each triangle not within a region will have a "0" added
to the end of its attribute list.

If the "a" switch is used without a number following, each triangle generated 
within the bounds of a region will be subject to that region's area
constraint.

If the "A" or "a" switches are not in effect, addRegion has the same effect as 
addPolygon.

=head1 METHODS TO ACCESS OUTPUT LISTS

The following methods retrieve the output lists after the triangulate method has
been invoked in void context.

Triangle's output data is not imported from C to Perl until one of these methods
is invoked, and then only what's needed to construct the list requested. So 
there may be a speed or memory advantage to accessing the output in this way - 
only what you need, when you need it.

The methods prefixed with "v" access the Voronoi diagram nodes and edges, if one
was generated.

=head2 nodes

Returns a reference to a list of nodes (vertices or points). 

    my $pointlist = $tri->nodes();    # retrieve nodes/vertices/points
    
The nodes in the list have this structure:

    [x, y, < zero or more attributes >, < boundary marker >]

=head2 elements

Returns a reference to a list of elements.

    $triangles  = $tri->elements(); # retrieve triangle list

The elements in the list have this structure:

    [[x0, y0], [x1, y1], [x2, y2],
     < another three vertices, if "o2" switch used >
     < zero or more attributes >
    ]

=head2 segments

Returns a reference to a list of segments.

    $segs  = $tri->segments(); # retrieve the PSLG segments

The segments in the list have this structure:

    [[x0, y0], [x1, y1], < boundary marker >]

=head2 edges

Returns a reference to a list of edges.

    $edges  = $tri->edges();    # retrieve all the triangle edges

The edges in the list have this structure:

    [[x0, y0], [x1, y1], < boundary marker >]

Note that the edge list is not produced by default. Request that it be generated
by invoking doEdges(1), or passing the 'e' switch to C<triangulate()>.

=head2 vnodes

Returns a reference to a list of nodes in the Voronoi diagram.

    $vpointlist = $tri->vnodes();   # retrieve Voronoi vertices

The Voronoi diagram nodes in the list have this structure:

    [x, y, < zero or more attributes >]

=head2 vedges

Returns a reference to a list of edges in the Voronoi diagram. Some of these
edges are actually rays.

    $vedges = $tri->vedges();   # retrieve Voronoi diagram edges and rays 

The Voronoi diagram edges in the list have this structure:

    [< vertex or vector >, < vertex or vector >, < ray flag >]

If the edge is a true edge, the ray flag will be 0.
If the edge is actually a ray, the ray flag will either be 1 or 2,
to indicate whether the the first, or second vertex should be interpreted as
a direction vector for the ray.

=head1 API STATUS

Currently Triangle's option strings are exposed to give more complete access to
its features. More of these options, and perhaps certain common combinations of 
them, will likely be wrapped in method-call getter-setters. I would prefer to 
preserve the ability to use the option strings directly, but it may be necessary
at some point to hide them completely behind a less idiosyncratic interface.


=head1 AUTHOR

Michael E. Sheldrake, C<< <sheldrake at cpan.org> >>

Triangle's author is Jonathan Richard Shewchuk


=head1 BUGS

Please report any bugs or feature requests to 

C<bug-math-geometry-delaunay at rt.cpan.org>

or through the web interface at 

L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Math-Geometry-Delaunay>

I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Math::Geometry::Delaunay


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Math-Geometry-Delaunay>

=item * Search CPAN

L<http://search.cpan.org/dist/Math-Geometry-Delaunay/>

=back


=head1 ACKNOWLEDGEMENTS

Thanks go to Far Leaves Tea in Berkeley for providing oolongs and refuge, 
and a place for paths to intersect.


=head1 LICENSE AND COPYRIGHT

Copyright 2012 Micheal E. Sheldrake.

This Perl binding to Triangle is free software; 
you can redistribute it and/or modify it under the terms of either: 
the GNU General Public License as published by the Free Software Foundation; 
or the Artistic License.

See http://dev.perl.org/licenses/ for more information.

=head2 Triangle license

B<Triangle> by Jonathan Richard Shewchuk, copyright 2005, includes the following
notice in the C source code. Please refer to the C source, included in with this
Perl module distribution, for the full notice.

    This program may be freely redistributed under the condition that the
    copyright notices (including this entire header and the copyright
    notice printed when the `-h' switch is selected) are not removed, and
    no compensation is received.  Private, research, and institutional
    use is free.  You may distribute modified versions of this code UNDER
    THE CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE TO IT IN THE
    SAME FILE REMAIN UNDER COPYRIGHT OF THE ORIGINAL AUTHOR, BOTH SOURCE
    AND OBJECT CODE ARE MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR
    NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution of this code as
    part of a commercial system is permissible ONLY BY DIRECT ARRANGEMENT
    WITH THE AUTHOR.  (If you are not directly supplying this code to a
    customer, and you are instead telling them how they can obtain it for
    free, then you are not required to make any arrangement with me.)

=cut

1;
