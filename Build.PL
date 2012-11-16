use strict;
use warnings;
use Module::Build;
use ExtUtils::CppGuess;

my $builder = Module::Build->new(
    module_name         => 'Math::Geometry::Delaunay',
    license             => 'perl',
    dist_author         => q{Michael E. Sheldrake <sheldrake@cpan.org>},
    configure_requires  => {'Module::Build'  => '0.38'},
    build_requires => {
        'Test::More' => 0,
        'ExtUtils::CBuilder' => 0,
        'ExtUtils::CppGuess' => 0
    },
    c_source            => ['src'],
    xs_files            => {'Delaunay.xs' => 'lib/Math/Geometry/Delaunay.xs'},
    
    # FLAGS FOR TRIANGLE
    # We want to include -DREDUCED, but it looks like that might take 
    # away the D option, even though the docs indicate otherwise. 
    # (Also, not sure the D option works, but need to test that more.)

    extra_compiler_flags => [qw(
                                -DREAL=double 
                                -DVOID=int 
                                -DTRILIBRARY 
                                -DANSI_DECLARATORS 
                                -DNO_TIMER
                               )
                            ],
    add_to_cleanup      => [ 'Math-Geometry-Delaunay-*' ],
);

# Triangle disables extended-precision doubles on x86 processors to make its
# exact arithmatic routines, which depend on the rounding error of normal-precision
# doubles, work properly. How to disable extended-precision doubles is 
# compiler-specific, and Triangle will do the right thing for x86 processors
# compiling with gcc under Linux or MSVC, if we pass either the LINUX
# or CPU86 compiler flag.

# However,
# we're also going to use xpfpa.h from 
# http://www.christian-seiler.de/projekte/fpmath/
# to try to set the FPU state, and more significantly, to restore the state
# to whateverever it was before each call to triangulate().
# The HAVE__CONTROLFP (Microsoft) and HAVE__FPU_SETCW (Linux) flags will give us
# the right macros for the platform.

my $compiler = ExtUtils::CppGuess->new();

if ($compiler->is_gcc) {
    if ($^O =~ /cygwin|MSWin32/i) { # MinGW or Cygwin - untested
        push @{$builder->extra_compiler_flags()}, "-DCPU86", "-DHAVE__CONTROLFP";
        }
    elsif ($^O =~ /linux/i) {
        push @{$builder->extra_compiler_flags()}, "-DLINUX", "-DHAVE__FPU_SETCW";
        }
    # else hope the system is in normal precision double mode by default
    }
elsif ($compiler->is_msvc) {
    push @{$builder->extra_compiler_flags()}, "-DCPU86", "-DHAVE__CONTROLFP";
    }

$builder->create_build_script();
