#!/bin/bash

# SeqSimulator Dependency Installation Script
# Supports Ubuntu/Debian and CentOS/RHEL systems

set -e

echo "========================================="
echo "SeqSimulator Dependency Installation"
echo "========================================="

# Detect OS
if [ -f /etc/os-release ]; then
    . /etc/os-release
    OS=$ID
else
    echo "Cannot detect OS. Please install dependencies manually."
    exit 1
fi

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to install on Ubuntu/Debian
install_debian() {
    echo "Detected Debian/Ubuntu system"
    
    # Update package list
    sudo apt-get update
    
    # Install basic dependencies
    sudo apt-get install -y \
        build-essential \
        cmake \
        git \
        wget \
        curl \
        python3 \
        python3-pip \
        python3-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        autoconf \
        automake
}

# Function to install on CentOS/RHEL
install_redhat() {
    echo "Detected RedHat/CentOS system"
    
    # Install development tools
    sudo yum groupinstall -y "Development Tools"
    
    # Install dependencies
    sudo yum install -y \
        cmake \
        git \
        wget \
        curl \
        python3 \
        python3-pip \
        python3-devel \
        zlib-devel \
        bzip2-devel \
        xz-devel \
        ncurses-devel \
        libcurl-devel \
        openssl-devel \
        autoconf \
        automake
}

# Install system packages
case "$OS" in
    ubuntu|debian)
        install_debian
        ;;
    centos|rhel|fedora)
        install_redhat
        ;;
    *)
        echo "Unsupported OS: $OS"
        echo "Please install dependencies manually"
        exit 1
        ;;
esac

# Create tools directory
TOOLS_DIR="$HOME/seqsim_tools"
mkdir -p "$TOOLS_DIR"
cd "$TOOLS_DIR"

echo ""
echo "Installing bioinformatics tools..."
echo ""

# Install SAMtools
if ! command_exists samtools; then
    echo "Installing SAMtools..."
    wget -q https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2
    tar -xjf samtools-1.17.tar.bz2
    cd samtools-1.17
    ./configure --prefix="$TOOLS_DIR"
    make -j$(nproc)
    make install
    cd ..
    rm -rf samtools-1.17*
else
    echo "SAMtools already installed"
fi

# Install BCFtools
if ! command_exists bcftools; then
    echo "Installing BCFtools..."
    wget -q https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2
    tar -xjf bcftools-1.17.tar.bz2
    cd bcftools-1.17
    ./configure --prefix="$TOOLS_DIR"
    make -j$(nproc)
    make install
    cd ..
    rm -rf bcftools-1.17*
else
    echo "BCFtools already installed"
fi

# Install BWA
if ! command_exists bwa; then
    echo "Installing BWA..."
    git clone https://github.com/lh3/bwa.git
    cd bwa
    make -j$(nproc)
    cp bwa "$TOOLS_DIR/bin/"
    cd ..
    rm -rf bwa
else
    echo "BWA already installed"
fi

# Install minimap2
if ! command_exists minimap2; then
    echo "Installing minimap2..."
    wget -q https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2
    tar -xjf minimap2-2.26_x64-linux.tar.bz2
    cp minimap2-2.26_x64-linux/minimap2 "$TOOLS_DIR/bin/"
    rm -rf minimap2-2.26*
else
    echo "minimap2 already installed"
fi

# Install ART
if ! command_exists art_illumina; then
    echo "Installing ART..."
    wget -q https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz
    tar -xzf artbinmountrainier2016.06.05linux64.tgz
    cp art_bin_MountRainier/art_illumina "$TOOLS_DIR/bin/"
    rm -rf art*
else
    echo "ART already installed"
fi

# Install DWGSIM
if ! command_exists dwgsim; then
    echo "Installing DWGSIM..."
    git clone --recursive https://github.com/nh13/DWGSIM.git
    cd DWGSIM
    make -j$(nproc)
    cp dwgsim "$TOOLS_DIR/bin/"
    cd ..
    rm -rf DWGSIM
else
    echo "DWGSIM already installed"
fi

# Install Python packages
echo ""
echo "Installing Python packages..."
pip3 install --user --upgrade pip
pip3 install --user \
    pyyaml \
    pysam \
    numpy \
    biopython \
    click \
    colorlog \
    tqdm \
    pandas \
    scipy \
    badread

# Install PBSIM3 (optional, larger download)
read -p "Install PBSIM3 for PacBio simulation? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    if ! command_exists pbsim3; then
        echo "Installing PBSIM3..."
        wget -q https://github.com/yukiteruono/pbsim3/releases/download/v3.0.0/pbsim3-3.0.0-Linux.tar.gz
        tar -xzf pbsim3-3.0.0-Linux.tar.gz
        cp pbsim3-3.0.0-Linux/pbsim3 "$TOOLS_DIR/bin/"
        rm -rf pbsim3*
    else
        echo "PBSIM3 already installed"
    fi
fi

# Add tools to PATH
echo ""
echo "========================================="
echo "Installation Complete!"
echo "========================================="
echo ""
echo "Add the following line to your ~/.bashrc or ~/.bash_profile:"
echo ""
echo "export PATH=\"$TOOLS_DIR/bin:\$PATH\""
echo ""
echo "Then run: source ~/.bashrc"
echo ""
echo "Or for immediate use, run:"
echo "export PATH=\"$TOOLS_DIR/bin:\$PATH\""
echo ""

# Verify installations
echo "Verifying installations..."
echo ""

# Check each tool
tools=("samtools" "bcftools" "bwa" "minimap2" "art_illumina" "dwgsim")
for tool in "${tools[@]}"; do
    if command_exists "$tool" || [ -f "$TOOLS_DIR/bin/$tool" ]; then
        echo "✓ $tool installed"
    else
        echo "✗ $tool not found"
    fi
done

# Check Python packages
echo ""
echo "Python packages:"
python3 -c "import yaml; print('✓ PyYAML installed')" 2>/dev/null || echo "✗ PyYAML not found"
python3 -c "import pysam; print('✓ pysam installed')" 2>/dev/null || echo "✗ pysam not found"
python3 -c "import numpy; print('✓ numpy installed')" 2>/dev/null || echo "✗ numpy not found"
python3 -c "import Bio; print('✓ biopython installed')" 2>/dev/null || echo "✗ biopython not found"
python3 -c "import badread; print('✓ badread installed')" 2>/dev/null || echo "✗ badread not found"

echo ""
echo "Installation script completed!"
