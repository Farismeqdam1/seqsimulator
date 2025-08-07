# SeqSimulator

A comprehensive, modular tool for simulating both short-read and long-read sequencing data with realistic coverage, error profiles, and reference-based variant generation.

## Features

- **Multi-platform Support**: Simulate Illumina short reads and ONT/PacBio long reads
- **Variant-aware Simulation**: Embed realistic genomic variants (SNPs, INDELs, SVs)
- **Multi-sample Generation**: Generate multiple samples with population-level variation
- **Flexible Configuration**: YAML-based configuration for easy customization
- **Complete Output**: FASTQ + BAM + VCF + documentation
- **Chromosome Subsampling**: Simulate specific chromosomes for faster testing
- **Docker Support**: Run in containerized environments

## Quick Start

### Installation

#### Method 1: Direct Installation
```bash
# Clone the repository
git clone https://github.com/Farismeqdam1/seqsimulator.git
cd seqsimulator

# Install Python dependencies
pip install -r requirements.txt

# Install required tools (Ubuntu/Debian)
sudo apt-get update
sudo apt-get install -y samtools bcftools bwa minimap2

# Download and install simulators
# ART (for short reads)
wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz
tar -xzf artbinmountrainier2016.06.05linux64.tgz
sudo cp art_bin_MountRainier/art_illumina /usr/local/bin/

# Badread (for long reads)
pip install badread
```

#### Method 2: Docker
```bash
# Build the Docker image
docker build -t seqsimulator .

# Run with Docker
docker run -v /path/to/data:/data seqsimulator -c /data/config.yaml
```

#### Method 3: Singularity
```bash
# Build Singularity image from Docker
singularity build seqsimulator.sif docker://seqsimulator

# Run with Singularity
singularity exec seqsimulator.sif python3 /app/seqsimulator.py -c config.yaml
```

### Basic Usage

1. **Generate example configuration:**
```bash
python seqsimulator.py --generate-config
```

2. **Edit the configuration file:**
```yaml
simulation_type: hybrid
reference_genome: /path/to/hg38.fa
output_directory: ./simulation_output
num_samples: 3

short_read_config:
  coverage: 30
  read_length: 150
  
long_read_config:
  coverage: 20
  mean_length: 10000
  platform: ont
```

3. **Run the simulation:**
```bash
python seqsimulator.py -c config.yaml -v
```

## Configuration Options

### Basic Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `simulation_type` | Type of simulation: `short_read`, `long_read`, or `hybrid` | Required |
| `reference_genome` | Path to reference FASTA file | Required |
| `output_directory` | Directory for output files | Required |
| `num_samples` | Number of samples to generate | 1 |
| `chromosomes` | List of chromosomes to simulate | All |
| `threads` | Number of threads to use | 4 |

### Short Read Configuration

```yaml
short_read_config:
  coverage: 30              # Target coverage (X)
  read_length: 150          # Read length in bp
  insert_size: 350          # Insert size for paired-end
  insert_sd: 50             # Insert size standard deviation
  
short_read_simulator: art   # Options: art, dwgsim
```

### Long Read Configuration

```yaml
long_read_config:
  coverage: 20              # Target coverage (X)
  mean_length: 10000        # Mean read length in bp
  error_rate: 0.1           # Base error rate
  platform: ont             # Platform: ont or pacbio
  
long_read_simulator: badread  # Options: badread, pbsim
```

### Variant Configuration

```yaml
variant_config:
  # Option 1: Use existing VCF
  population_vcf: /path/to/1000genomes.vcf.gz
  
  # Option 2: Generate random variants
  snp_rate: 0.001           # SNPs per base
  indel_rate: 0.0001        # INDELs per base
  sv_rate: 0.00001          # SVs per base
```

## Output Structure

```
simulation_output/
├── fastq/                  # Sequencing reads
│   ├── sample_001_R1.fastq.gz
│   ├── sample_001_R2.fastq.gz
│   └── sample_001_long.fastq.gz
├── bam/                    # Aligned reads
│   ├── sample_001.sorted.bam
│   └── sample_001.sorted.bam.bai
├── vcf/                    # Variant files
│   └── sample_001.vcf.gz
├── logs/                   # Simulation logs
├── reports/                # Summary reports
│   └── summary.json
└── README.md               # Simulation details
```

## Examples

### Example 1: Simulate Illumina Whole Genome Sequencing
```yaml
simulation_type: short_read
reference_genome: /ref/hg38.fa
output_directory: ./wgs_simulation
num_samples: 5

short_read_config:
  coverage: 30
  read_length: 150
  insert_size: 350

variant_config:
  snp_rate: 0.001
  indel_rate: 0.0001
```

### Example 2: Simulate ONT Long Reads for Chr22
```yaml
simulation_type: long_read
reference_genome: /ref/hg38.fa
output_directory: ./ont_chr22
chromosomes: [chr22]

long_read_config:
  coverage: 40
  mean_length: 15000
  platform: ont
  error_rate: 0.05
```

### Example 3: Hybrid Sequencing with Population Variants
```yaml
simulation_type: hybrid
reference_genome: /ref/hg38.fa
output_directory: ./hybrid_pop
num_samples: 10

short_read_config:
  coverage: 30
  read_length: 150

long_read_config:
  coverage: 15
  mean_length: 10000

variant_config:
  population_vcf: /data/1000G_phase3.vcf.gz
```

## Integration with Downstream Tools

### GATK Variant Calling
```bash
# Index reference
samtools faidx reference.fa
gatk CreateSequenceDictionary -R reference.fa

# Call variants
gatk HaplotypeCaller \
    -R reference.fa \
    -I simulation_output/bam/sample_001.sorted.bam \
    -O calls.vcf
```

### Exomiser Analysis
```bash
# Prepare PED file
echo -e "FAM001\tsample_001\t0\t0\t1\t2" > family.ped

# Run Exomiser
java -jar exomiser.jar \
    --analysis-config-file analysis.yml \
    --vcf simulation_output/vcf/sample_001.vcf.gz \
    --ped family.ped
```

### vg Graph Construction
```bash
# Build variation graph
vg construct \
    -r reference.fa \
    -v simulation_output/vcf/sample_001.vcf.gz \
    > graph.vg

# Index the graph
vg index -x graph.xg graph.vg
```

## Advanced Features

### Custom Error Profiles
Create custom error profiles for specific sequencing platforms:

```python
# custom_profile.py
from seqsimulator import ErrorProfile

profile = ErrorProfile()
profile.substitution_rate = 0.001
profile.insertion_rate = 0.0005
profile.deletion_rate = 0.0005
profile.save("custom_ont_profile.json")
```

### Pedigree-based Simulation
Simulate families with inheritance patterns:

```yaml
pedigree_config:
  ped_file: family.ped
  inheritance_mode: mendelian
  de_novo_rate: 1e-8
```

### Quality Score Calibration
Calibrate quality scores based on real data:

```bash
python scripts/calibrate_quality.py \
    --real-fastq real_data.fastq.gz \
    --output quality_model.json
```

## Performance Optimization

### Memory Usage
- For large genomes, use chromosome subsampling
- Adjust coverage to reduce memory footprint
- Use `--low-memory` flag for resource-constrained systems

### Speed Optimization
- Increase thread count: `threads: 16`
- Use faster simulators (ART over DWGSIM for short reads)
- Enable parallel sample generation

### Storage Optimization
- Compress outputs: All outputs are gzipped by default
- Remove intermediate files: `keep_intermediate: false`
- Use symbolic links for reference genome

## Troubleshooting

### Common Issues

**Issue**: "Reference genome not indexed"
```bash
# Solution: Index the reference
samtools faidx reference.fa
bwa index reference.fa
```

**Issue**: "Tool not found in PATH"
```bash
# Solution: Install missing tool or use Docker
sudo apt-get install samtools
# OR
docker run -v $(pwd):/data seqsimulator -c /data/config.yaml
```

**Issue**: "Out of memory error"
```bash
# Solution: Reduce coverage or use chromosome subsampling
chromosomes: [chr21, chr22]  # Simulate only small chromosomes
```

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Development Setup
```bash
# Clone repository
git clone https://github.com/yourusername/seqsimulator.git
cd seqsimulator

# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install in development mode
pip install -e .

# Run tests
pytest tests/
```

## Citation

If you use SeqSimulator in your research, please cite:

```bibtex
@software{seqsimulator2024,
  title = {SeqSimulator: A modular tool for sequencing data simulation},
  author = {SeqSim Development Team},
  year = {2024},
  url = {https://github.com/Farismeqdam1/seqsimulator}
}
```

