#!/usr/bin/env bash

# Unified container build script
# Usage: 
#   ./build_containers.sh                    # Build both Docker and Singularity (default)
#   ./build_containers.sh --docker           # Build only Docker images
#   ./build_containers.sh --singularity      # Build only Singularity images
#   ./build_containers.sh --all              # Build both (explicit)
#   ./build_containers.sh --github_action      # Build both for GitHub Actions
#   ./build_containers.sh --docker github_action    # Build Docker for GitHub Actions
#   ./build_containers.sh --singularity github_action # Build Singularity for GitHub Actions

set -e

# Parse arguments
build_docker=false
build_singularity=false
github_action_mode=""
auto_detect=false

# If no arguments provided, auto-detect available tools
if [ $# -eq 0 ]; then
    auto_detect=true
fi

# Parse all arguments
while [ $# -gt 0 ]; do
    case "$1" in
        --docker)
            build_docker=true
            ;;
        --singularity|--apptainer)
            build_singularity=true
            ;;
        --all)
            build_docker=true
            build_singularity=true
            ;;
        --github_action)
            github_action_mode="github_action"
            # If --github_action is used without --docker or --singularity, enable auto-detect
            if [ "$build_docker" = false ] && [ "$build_singularity" = false ]; then
                auto_detect=true
            fi
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS] [github_action]"
            echo ""
            echo "Options:"
            echo "  --docker          Build Docker images only"
            echo "  --singularity     Build Singularity/Apptainer images only"
            echo "  --apptainer       Alias for --singularity"
            echo "  --all             Build both Docker and Singularity images"
            echo "  (no options)      Build both by default"
            echo ""
            echo "Additional arguments:"
            echo "  github_action     Enable GitHub Actions mode (saves to archive)"
            echo ""
            echo "Examples:"
            echo "  $0                              # Build both (default)"
            echo "  $0 --docker                     # Build only Docker"
            echo "  $0 --singularity github_action  # Build Singularity for CI"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
    shift
done

# Auto-detect available tools if no specific option was provided
if [ "$auto_detect" = true ]; then
    echo "Auto-detecting available container tools..."
    
    if command -v docker &> /dev/null; then
        echo "  ✓ Docker found"
        build_docker=true
    else
        echo "  ✗ Docker not found"
    fi
    
    if command -v singularity &> /dev/null; then
        echo "  ✓ Singularity/Apptainer found"
        build_singularity=true
    else
        echo "  ✗ Singularity/Apptainer not found"
    fi
    
    # Check if at least one tool is available
    if [ "$build_docker" = false ] && [ "$build_singularity" = false ]; then
        echo ""
        echo "❌ Error: Neither Docker nor Singularity/Apptainer found on this system."
        echo "Please install at least one container tool to proceed."
        exit 1
    fi
    
    echo ""
fi

# Save original working directory
wd=$(pwd)

# Get architecture for Docker builds
arch=$(uname -m)

# Build Docker images if requested
if [ "$build_docker" = true ]; then
    echo "════════════════════════════════════════════════════════════════"
    echo "  Building Docker images..."
    echo "════════════════════════════════════════════════════════════════"
    
    # List of image names
    image_list=( )
    
    for dir in docker/*; do 
        cd "${dir}"
        imgname=$(echo $dir | rev | cut -d/ -f1 | rev)
        image_list+=(${imgname})

        echo ██████████████████▓▒░   Building ${imgname}   ░▒▓██████████████████
        
        # Set architecture to docker buildx
        docker_arch_option=""
        
        # Reditools2 and Pluviometer do not compile on arm64, force using amd64 compilation
        if [[ "$arch" == arm* || "$arch" == "aarch64" ]]; then
            if [[ $dir =~ "reditools2" ]] || [[ $dir =~ "pluviometer" ]]; then
                echo "$imgname does not compile on arm64, force using amd64 compilation"
                docker_arch_option=" --platform linux/amd64"
            fi
        fi

        docker build ${docker_arch_option} -t ${imgname} .
        
        # Back to the original working directory
        cd "$wd"
    done

    if [[ ${github_action_mode} == 'github_action' ]]; then
        echo "Saving docker images to cache..."
        docker save ${image_list[@]} -o docker-images.tar
        echo Archive size: $(stat --printf="%s" docker-images.tar)
    fi
    
    echo ""
fi

# Build Singularity images if requested
if [ "$build_singularity" = true ]; then
    echo "════════════════════════════════════════════════════════════════"
    echo "  Building Singularity/Apptainer images..."
    echo "════════════════════════════════════════════════════════════════"
    
    # Create output directory for .sif images
    mkdir -p sif_images

    # List to store built images
    sif_list=( )

    # Loop through all .def files in singularity/ subdirectories
    for def_file in singularity/*/*.def; do
        # Extract image name from directory
        imgname=$(basename "$(dirname "$def_file")")
        sif_output="sif_images/${imgname}.sif"
        sif_list+=("$sif_output")

        echo ██████████████████▓▒░   Building ${imgname}   ░▒▓██████████████████
        
        # Build with fakeroot if available
        if command -v singularity &> /dev/null; then
            singularity build --fakeroot "$sif_output" "$def_file"
        else
            echo "❌ Singularity not found! Aborting."
            exit 1
        fi
    done

    # If used in GitHub Actions, create an archive
    if [[ ${github_action_mode} == 'github_action' ]]; then
        echo "Saving SIF images to cache archive..."
        tar -cf singularity-images.tar "${sif_list[@]}"
        echo Archive size: $(stat --printf="%s" singularity-images.tar)
    fi
    
    echo ""
fi

echo "════════════════════════════════════════════════════════════════"
echo "  ✓ Build complete!"
echo "════════════════════════════════════════════════════════════════"
