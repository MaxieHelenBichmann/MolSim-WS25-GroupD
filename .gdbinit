# ==============================================================================
# GDB Configuration for MolSim Project
# ==============================================================================
# This file is automatically loaded when GDB starts in this directory.
# Use: gdb ./build/project/MolSim
#
# IMPORTANT: If GDB gives a security warning, add this to ~/.gdbinit:
#   add-auto-load-safe-path /home/henry/MolSim/MolSim-WS25-GroupD/.gdbinit
# ==============================================================================

# ------------------------------------------------------------------------------
# General Settings
# ------------------------------------------------------------------------------

# Enable pretty printing for STL containers (vector, array, etc.)
set print pretty on
set print array on
set print array-indexes on
set print elements 200

# Show full object details without truncation
set print object on
set print static-members on
set print vtbl on

# Demangle C++ symbols for readability
set print demangle on
set demangle-style gnu-v3

# Show more context when stopped
set listsize 20

# Enable history and save command history
set history save on
set history size 10000
set history filename ~/.gdb_history

# Auto-confirm dangerous operations (comment out if you prefer manual confirmation)
# set confirm off

# Enable TUI mode helpers
set tui border-kind ascii
set tui active-border-mode bold-standout

# Better thread display
set print thread-events on

# ------------------------------------------------------------------------------
# C++20 and Template Debugging
# ------------------------------------------------------------------------------

# Skip standard library stepping (speeds up debugging)
skip -gfi /usr/include/c++/*
skip -gfi /usr/include/*/c++/*
skip -gfi /usr/lib/gcc/*

# ------------------------------------------------------------------------------
# MolSim-Specific Pretty Printers
# ------------------------------------------------------------------------------

# Custom command to print a Vector<T, N> (mol_sim::R3 or mol_sim::Vector)
define pvec
    if $argc == 0
        help pvec
    else
        printf "Vector<%s, %zu>: [", $arg0.data_._M_elems[0]@1, sizeof($arg0.data_._M_elems)/sizeof($arg0.data_._M_elems[0])
        set $i = 0
        set $size = sizeof($arg0.data_._M_elems)/sizeof($arg0.data_._M_elems[0])
        while $i < $size
            if $i > 0
                printf ", "
            end
            printf "%.6f", $arg0.data_._M_elems[$i]
            set $i = $i + 1
        end
        printf "]\n"
    end
end
document pvec
Print a mol_sim::Vector<T, N> in a readable format.
Usage: pvec <vector_variable>
Example: pvec particle.x
end

# Custom command to print a Particle with all its attributes
define pparticle
    if $argc == 0
        help pparticle
    else
        printf "Particle @ %p:\n", &$arg0
        printf "  Position (x):    "
        pvec $arg0.x
        printf "  Old Position:    "
        pvec $arg0.old_x
        printf "  Velocity (v):    "
        pvec $arg0.v
        printf "  Force (f):       "
        pvec $arg0.f
        printf "  Old Force:       "
        pvec $arg0.old_f
        printf "  Mass (m):        %.6f\n", $arg0.m
        printf "  Epsilon:         %.6f\n", $arg0.epsilon
        printf "  Sigma:           %.6f\n", $arg0.sigma
        printf "  Type:            %d\n", $arg0.type
    end
end
document pparticle
Print a mol_sim::Particle with all its attributes in a formatted way.
Usage: pparticle <particle_variable>
Example: pparticle p1
         pparticle particles[0]
end

# Print a particle container summary
define pcontainer
    if $argc == 0
        help pcontainer
    else
        printf "ParticleContainer @ %p:\n", &$arg0
        printf "  Size: %zu particles\n", $arg0.size()
        if $argc > 1
            printf "  First %d particles:\n", $arg1
            set $i = 0
            while $i < $arg1 && $i < $arg0.size()
                printf "  [%d] ", $i
                pparticle $arg0[$i]
                printf "\n"
                set $i = $i + 1
            end
        end
    end
end
document pcontainer
Print a ParticleContainer summary and optionally the first N particles.
Usage: pcontainer <container> [count]
Example: pcontainer particles 5
end

# Print R3 vector (common alias)
define pr3
    if $argc == 0
        help pr3
    else
        pvec $arg0
    end
end
document pr3
Shorthand for printing a mol_sim::R3 vector (alias for pvec).
Usage: pr3 <r3_variable>
Example: pr3 particle.getX()
end

# ------------------------------------------------------------------------------
# Breakpoint Shortcuts for Common MolSim Functions
# ------------------------------------------------------------------------------

# Break at simulation start
define bsim
    break mol_sim::Simulation::run
end
document bsim
Set breakpoint at the start of Simulation::run().
Usage: bsim
end

# Break at force calculation
define bforce
    break mol_sim::Simulation::calculateF
end
document bforce
Set breakpoint at force calculation in Simulation::calculateF().
Usage: bforce
end

# Break at position update
define bpos
    break mol_sim::Simulation::calculateX
end
document bpos
Set breakpoint at position calculation in Simulation::calculateX().
Usage: bpos
end

# Break at velocity update
define bvel
    break mol_sim::Simulation::calculateV
end
document bvel
Set breakpoint at velocity calculation in Simulation::calculateV().
Usage: bvel
end

# Break at boundary application
define bboundary
    break mol_sim::Simulation::applyBoundaries
end
document bboundary
Set breakpoint at boundary condition application in Simulation::applyBoundaries().
Usage: bboundary
end

# Break at main
define bmain
    break main
end
document bmain
Set breakpoint at main function.
Usage: bmain
end

# ------------------------------------------------------------------------------
# Watchpoint Helpers
# ------------------------------------------------------------------------------

# Watch a particle's position
define watchpos
    if $argc == 0
        help watchpos
    else
        watch $arg0.x.data_._M_elems[0]
        watch $arg0.x.data_._M_elems[1]
        watch $arg0.x.data_._M_elems[2]
    end
end
document watchpos
Set watchpoints on all components of a particle's position vector.
Usage: watchpos <particle>
Example: watchpos particles[0]
end

# Watch a particle's velocity
define watchvel
    if $argc == 0
        help watchvel
    else
        watch $arg0.v.data_._M_elems[0]
        watch $arg0.v.data_._M_elems[1]
        watch $arg0.v.data_._M_elems[2]
    end
end
document watchvel
Set watchpoints on all components of a particle's velocity vector.
Usage: watchvel <particle>
Example: watchvel particles[0]
end

# ------------------------------------------------------------------------------
# Analysis Commands
# ------------------------------------------------------------------------------

# Calculate total kinetic energy of all particles in a container
define kenergy
    if $argc == 0
        help kenergy
    else
        set $ke = 0.0
        set $i = 0
        while $i < $arg0.size()
            set $p = $arg0[$i]
            set $v2 = $p.v.data_._M_elems[0]*$p.v.data_._M_elems[0] + \
                      $p.v.data_._M_elems[1]*$p.v.data_._M_elems[1] + \
                      $p.v.data_._M_elems[2]*$p.v.data_._M_elems[2]
            set $ke = $ke + 0.5 * $p.m * $v2
            set $i = $i + 1
        end
        printf "Total kinetic energy: %.6e\n", $ke
    end
end
document kenergy
Calculate and print the total kinetic energy of all particles in a container.
Usage: kenergy <container>
Example: kenergy particles
end

# Print center of mass
define pcom
    if $argc == 0
        help pcom
    else
        set $com_x = 0.0
        set $com_y = 0.0
        set $com_z = 0.0
        set $total_m = 0.0
        set $i = 0
        while $i < $arg0.size()
            set $p = $arg0[$i]
            set $com_x = $com_x + $p.m * $p.x.data_._M_elems[0]
            set $com_y = $com_y + $p.m * $p.x.data_._M_elems[1]
            set $com_z = $com_z + $p.m * $p.x.data_._M_elems[2]
            set $total_m = $total_m + $p.m
            set $i = $i + 1
        end
        printf "Center of mass: [%.6f, %.6f, %.6f]\n", $com_x/$total_m, $com_y/$total_m, $com_z/$total_m
    end
end
document pcom
Calculate and print the center of mass of all particles in a container.
Usage: pcom <container>
Example: pcom particles
end

# Calculate current temperature
define ptemp
    if $argc < 2
        help ptemp
    else
        set $ke = 0.0
        set $i = 0
        while $i < $arg0.size()
            set $p = $arg0[$i]
            set $v2 = $p.v.data_._M_elems[0]*$p.v.data_._M_elems[0] + \
                      $p.v.data_._M_elems[1]*$p.v.data_._M_elems[1] + \
                      $p.v.data_._M_elems[2]*$p.v.data_._M_elems[2]
            set $ke = $ke + 0.5 * $p.m * $v2
            set $i = $i + 1
        end
        set $temp = (2.0 * $ke) / ($arg0.size() * $arg1)
        printf "Current temperature: %.6e (KE=%.6e, N=%zu, dim=%d)\n", $temp, $ke, $arg0.size(), $arg1
    end
end
document ptemp
Calculate and print the current temperature of the system.
Temperature is calculated as: T = (2 * kinetic_energy) / (N * dimensions)
Usage: ptemp <container> <dimensions>
Example: ptemp particles 3
         ptemp particles 2
end



# ------------------------------------------------------------------------------
# Debugging Workflow Shortcuts
# ------------------------------------------------------------------------------

# Run until next iteration of simulation loop
define nextiter
    break mol_sim::Simulation::calculateX
    continue
    delete breakpoints $bpnum
end
document nextiter
Continue execution until the next iteration of the simulation loop.
Usage: nextiter
end

# Print current simulation state
define psimstate
    if $argc == 0
        help psimstate
    else
        printf "=== Simulation State ===\n"
        printf "Current time:     %.6f\n", $arg0.current_time
        printf "Iteration:        %d\n", $arg0.iteration
        printf "Delta t:          %.6f\n", $arg0.delta_t
        printf "Particles:        %zu\n", $arg0.particles.size()
        printf "Total energy:     %.6e\n", $arg0.total_energy
        printf "=======================\n"
    end
end
document psimstate
Print the current state of a Simulation object.
Usage: psimstate <simulation>
end

# ------------------------------------------------------------------------------
# TUI Layout Presets
# ------------------------------------------------------------------------------

define layout-src
    layout src
    focus cmd
end
document layout-src
Switch to TUI source layout with command focus.
Usage: layout-src
end

define layout-asm
    layout asm
    focus cmd
end
document layout-asm
Switch to TUI assembly layout with command focus.
Usage: layout-asm
end

define layout-split
    layout split
    focus cmd
end
document layout-split
Switch to TUI split (source + assembly) layout with command focus.
Usage: layout-split
end

# ------------------------------------------------------------------------------
# Catchpoints for Exception Debugging
# ------------------------------------------------------------------------------

# Catch all C++ exceptions
define catchall
    catch throw
    catch catch
end
document catchall
Set catchpoints for all C++ exceptions (both throw and catch).
Usage: catchall
end

# ------------------------------------------------------------------------------
# Project-Specific Initialization
# ------------------------------------------------------------------------------

# Set breakpoint at common exception classes
# catch throw mol_sim::MolSimException
# catch throw mol_sim::SimulationException
# catch throw std::runtime_error

# Print welcome message
printf "\n"
printf "====================================================================\n"
printf "  MolSim GDB Configuration Loaded\n"
printf "====================================================================\n"
printf "Custom commands available:\n"
printf "  pvec, pr3         - Print Vector<T,N> or R3\n"
printf "  pparticle         - Print Particle with all attributes\n"
printf "  pcontainer        - Print ParticleContainer summary\n"
printf "  bsim, bforce, ... - Set breakpoints at key simulation functions\n"
printf "  watchpos, watchvel - Watch particle position/velocity changes\n"
printf "  kenergy, pcom, ptemp - Calculate energy, center of mass, temperature\n"
printf "  nextiter          - Continue to next simulation iteration\n"
printf "  catchall          - Catch all C++ exceptions\n"
printf "\n"
printf "Type 'help <command>' for detailed usage information.\n"
printf "====================================================================\n"
printf "\n"

# Auto-load symbol file if it exists (adjust path as needed)
# symbol-file build/project/MolSim

# Set initial breakpoint at main (comment out if not desired)
# break main
