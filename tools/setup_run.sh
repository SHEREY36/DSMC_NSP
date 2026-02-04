#!/bin/bash
# setup_run.sh - Create a new DSMC run directory with proper structure
# Usage: ./setup_run.sh <run_name>

set -e

if [ $# -ne 1 ]; then
    echo "Usage: $0 <run_name>"
    echo "Example: $0 my_shear_case"
    exit 1
fi

RUN_NAME=$1
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
RUN_DIR="$PROJECT_ROOT/runs/$RUN_NAME"

# Check if run directory already exists
if [ -d "$RUN_DIR" ]; then
    echo "Error: Run directory already exists: $RUN_DIR"
    exit 1
fi

# Create run directory
echo "Creating run directory: $RUN_DIR"
mkdir -p "$RUN_DIR"

# Copy template input file
if [ -f "$PROJECT_ROOT/runs/template/system_input.dat" ]; then
    cp "$PROJECT_ROOT/runs/template/system_input.dat" "$RUN_DIR/"
    echo "Copied template system_input.dat"

    # Also copy the commented reference file
    if [ -f "$PROJECT_ROOT/runs/template/system_input_commented.dat" ]; then
        cp "$PROJECT_ROOT/runs/template/system_input_commented.dat" "$RUN_DIR/"
        echo "Copied system_input_commented.dat (reference with explanations)"
    fi
else
    echo "Warning: Template system_input.dat not found, creating basic template"
    cat > "$RUN_DIR/system_input.dat" << 'EOF'
0.01 0.01 0.01
10 10 10
10 10 10
1000
1.0
0.001 2.0
2500.0
0.9
1.0e-5 0.1
T T T
0.0 0.0 0.0
0.0 0.0 0.0
0.0 0.0 0.0
F
EOF
fi

# Create symlink to executable
cd "$RUN_DIR"
ln -sf ../../bin/DSMC ./DSMC
echo "Created symlink to executable"

# Create symlink to collision_models
ln -sf ../../collision_models ./collision_models
echo "Created symlink to collision_models/"

# Success message
echo ""
echo "========================================"
echo "Run directory created successfully!"
echo "========================================"
echo "Location: $RUN_DIR"
echo ""
echo "Next steps:"
echo "1. cd $RUN_DIR"
echo "2. Edit system_input.dat with your simulation parameters"
echo "3. Copy or generate particle_input.dat"
echo "4. Run: ./DSMC"
echo ""
echo "The simulation will automatically load the correct GMM file"
echo "based on the AR value in system_input.dat"
echo ""
