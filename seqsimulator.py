#!/usr/bin/env python3
"""
SeqSimulator - A modular tool for simulating short-read and long-read sequencing data
Author: SeqSim Development Team
Version: 1.0.0
"""

import argparse
import yaml
import os
import sys
import subprocess
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import json
import random
import shutil
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('SeqSimulator')


class SeqSimulator:
    """Main class for sequencing simulation"""
    
    def __init__(self, config_path: str):
        """Initialize the simulator with configuration"""
        self.config = self.load_config(config_path)
        self.validate_config()
        self.setup_output_directory()
        self.tool_paths = self.check_dependencies()
        
    def load_config(self, config_path: str) -> Dict:
        """Load YAML configuration file"""
        try:
            with open(config_path, 'r') as f:
                config = yaml.safe_load(f)
            logger.info(f"Configuration loaded from {config_path}")
            return config
        except Exception as e:
            logger.error(f"Failed to load configuration: {e}")
            sys.exit(1)
            
    def validate_config(self):
        """Validate configuration parameters"""
        required_fields = ['simulation_type', 'reference_genome', 'output_directory']
        for field in required_fields:
            if field not in self.config:
                logger.error(f"Missing required field: {field}")
                sys.exit(1)
                
        # Validate simulation type
        valid_types = ['short_read', 'long_read', 'hybrid']
        if self.config['simulation_type'] not in valid_types:
            logger.error(f"Invalid simulation type. Must be one of: {valid_types}")
            sys.exit(1)
            
        # Check if reference genome exists
        if not os.path.exists(self.config['reference_genome']):
            logger.error(f"Reference genome not found: {self.config['reference_genome']}")
            sys.exit(1)
            
    def setup_output_directory(self):
        """Create output directory structure"""
        output_dir = Path(self.config['output_directory'])
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create subdirectories
        subdirs = ['fastq', 'bam', 'vcf', 'logs', 'reports']
        for subdir in subdirs:
            (output_dir / subdir).mkdir(exist_ok=True)
            
        logger.info(f"Output directory structure created at {output_dir}")
        
    def check_dependencies(self) -> Dict[str, str]:
        """Check for required tools and return their paths"""
        tools = {}
        required_tools = []
        
        if self.config['simulation_type'] in ['short_read', 'hybrid']:
            if self.config.get('short_read_simulator', 'art') == 'art':
                required_tools.append('art_illumina')
            elif self.config.get('short_read_simulator') == 'dwgsim':
                required_tools.append('dwgsim')
                
        if self.config['simulation_type'] in ['long_read', 'hybrid']:
            if self.config.get('long_read_simulator', 'pbsim') == 'pbsim':
                required_tools.append('pbsim3')
            elif self.config.get('long_read_simulator') == 'badread':
                required_tools.append('badread')
                
        # Common tools
        required_tools.extend(['samtools', 'bcftools', 'bgzip', 'tabix'])
        
        for tool in required_tools:
            tool_path = shutil.which(tool)
            if tool_path:
                tools[tool] = tool_path
                logger.info(f"Found {tool}: {tool_path}")
            else:
                logger.warning(f"Tool {tool} not found in PATH")
                
        return tools
        
    def simulate_variants(self, sample_id: str) -> str:
        """Generate variants for a sample"""
        logger.info(f"Generating variants for sample {sample_id}")
        
        vcf_dir = Path(self.config['output_directory']) / 'vcf'
        output_vcf = vcf_dir / f"{sample_id}.vcf.gz"
        
        if 'variant_config' in self.config:
            var_config = self.config['variant_config']
            
            # Use bcftools to generate random variants
            cmd = [
                'bcftools', 'view',
                '-O', 'z',
                '-o', str(output_vcf)
            ]
            
            if 'population_vcf' in var_config:
                # Use existing population VCF
                cmd.extend(['-s', sample_id, var_config['population_vcf']])
            else:
                # Generate random variants
                self.generate_random_variants(sample_id, output_vcf, var_config)
                return str(output_vcf)
                
            try:
                subprocess.run(cmd, check=True)
                subprocess.run(['tabix', '-p', 'vcf', str(output_vcf)], check=True)
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to generate variants: {e}")
                
        return str(output_vcf)
        
    def generate_random_variants(self, sample_id: str, output_vcf: Path, var_config: Dict):
        """Generate random variants based on configuration"""
        import pysam
        
        ref = pysam.FastaFile(self.config['reference_genome'])
        
        # Parameters
        snp_rate = var_config.get('snp_rate', 0.001)
        indel_rate = var_config.get('indel_rate', 0.0001)
        sv_rate = var_config.get('sv_rate', 0.00001)
        
        vcf_header = self.create_vcf_header(sample_id, ref.references)
        
        variants = []
        
        # Select chromosomes to process
        chromosomes = self.config.get('chromosomes', ref.references)
        
        for chrom in chromosomes:
            if chrom not in ref.references:
                continue
                
            chrom_seq = ref.fetch(chrom)
            chrom_len = len(chrom_seq)
            
            # Generate SNPs
            num_snps = int(chrom_len * snp_rate)
            for _ in range(num_snps):
                pos = random.randint(1, chrom_len - 1)
                ref_base = chrom_seq[pos - 1].upper()
                if ref_base in 'ACGT':
                    alt_base = random.choice([b for b in 'ACGT' if b != ref_base])
                    qual = random.uniform(20, 60)
                    gt = random.choice(['0/1', '1/1'])
                    variants.append((chrom, pos, ref_base, alt_base, qual, gt, 'SNP'))
                    
            # Generate INDELs
            num_indels = int(chrom_len * indel_rate)
            for _ in range(num_indels):
                pos = random.randint(1, chrom_len - 10)
                indel_len = random.randint(1, 10)
                
                if random.random() < 0.5:  # Deletion
                    ref_seq = chrom_seq[pos - 1:pos + indel_len].upper()
                    alt_seq = ref_seq[0]
                else:  # Insertion
                    ref_seq = chrom_seq[pos - 1].upper()
                    alt_seq = ref_seq + ''.join(random.choices('ACGT', k=indel_len))
                    
                qual = random.uniform(20, 50)
                gt = random.choice(['0/1', '1/1'])
                variants.append((chrom, pos, ref_seq, alt_seq, qual, gt, 'INDEL'))
                
        # Sort variants
        variants.sort(key=lambda x: (x[0], x[1]))
        
        # Write VCF
        with open(output_vcf, 'w') as f:
            f.write(vcf_header)
            for var in variants:
                chrom, pos, ref, alt, qual, gt, var_type = var
                dp = random.randint(10, 100)
                info = f"TYPE={var_type};DP={dp}"
                format_field = "GT:DP"
                sample_field = f"{gt}:{dp}"
                f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual:.2f}\tPASS\t{info}\t{format_field}\t{sample_field}\n")
                
        # Compress and index
        subprocess.run(['bgzip', '-f', str(output_vcf)], check=True)
        subprocess.run(['tabix', '-p', 'vcf', f"{output_vcf}.gz"], check=True)
        
        ref.close()
        
    def create_vcf_header(self, sample_id: str, contigs: List[str]) -> str:
        """Create VCF header"""
        header = "##fileformat=VCFv4.3\n"
        header += f"##fileDate={datetime.now().strftime('%Y%m%d')}\n"
        header += f"##source=SeqSimulator_v1.0.0\n"
        header += f"##reference={self.config['reference_genome']}\n"
        
        for contig in contigs:
            header += f"##contig=<ID={contig}>\n"
            
        header += '##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type">\n'
        header += '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">\n'
        header += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        header += '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Sample depth">\n'
        header += f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}\n"
        
        return header
        
    def simulate_short_reads(self, sample_id: str, vcf_file: Optional[str] = None):
        """Simulate short-read sequencing data"""
        logger.info(f"Simulating short reads for sample {sample_id}")
        
        output_dir = Path(self.config['output_directory'])
        fastq_dir = output_dir / 'fastq'
        
        # Get parameters
        sr_config = self.config.get('short_read_config', {})
        coverage = sr_config.get('coverage', 30)
        read_length = sr_config.get('read_length', 150)
        insert_size = sr_config.get('insert_size', 350)
        insert_sd = sr_config.get('insert_sd', 50)
        
        simulator = self.config.get('short_read_simulator', 'art')
        
        if simulator == 'art' and 'art_illumina' in self.tool_paths:
            # Prepare reference (with variants if provided)
            ref_to_use = self.config['reference_genome']
            if vcf_file:
                ref_to_use = self.apply_variants_to_reference(sample_id, vcf_file)
                
            cmd = [
                self.tool_paths['art_illumina'],
                '-ss', 'HS25',  # Illumina HiSeq 2500
                '-i', ref_to_use,
                '-l', str(read_length),
                '-f', str(coverage),
                '-m', str(insert_size),
                '-s', str(insert_sd),
                '-p',  # Paired-end
                '-o', str(fastq_dir / f"{sample_id}_R")
            ]
            
            if self.config.get('chromosomes'):
                # ART doesn't support chromosome selection directly
                # Would need to extract chromosomes first
                pass
                
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                logger.info(f"Short reads generated for {sample_id}")
                
                # Rename output files
                os.rename(fastq_dir / f"{sample_id}_R1.fq", 
                         fastq_dir / f"{sample_id}_R1.fastq")
                os.rename(fastq_dir / f"{sample_id}_R2.fq", 
                         fastq_dir / f"{sample_id}_R2.fastq")
                
                # Compress FASTQ files
                subprocess.run(['gzip', str(fastq_dir / f"{sample_id}_R1.fastq")], check=True)
                subprocess.run(['gzip', str(fastq_dir / f"{sample_id}_R2.fastq")], check=True)
                
                # Generate BAM file
                self.generate_bam(sample_id, 'short')
                
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to generate short reads: {e}")
                
        elif simulator == 'dwgsim' and 'dwgsim' in self.tool_paths:
            self.simulate_with_dwgsim(sample_id, vcf_file)
        else:
            logger.error(f"Simulator {simulator} not available")
            
    def simulate_long_reads(self, sample_id: str, vcf_file: Optional[str] = None):
        """Simulate long-read sequencing data"""
        logger.info(f"Simulating long reads for sample {sample_id}")
        
        output_dir = Path(self.config['output_directory'])
        fastq_dir = output_dir / 'fastq'
        
        # Get parameters
        lr_config = self.config.get('long_read_config', {})
        coverage = lr_config.get('coverage', 20)
        mean_length = lr_config.get('mean_length', 10000)
        error_rate = lr_config.get('error_rate', 0.1)
        platform = lr_config.get('platform', 'ont')
        
        simulator = self.config.get('long_read_simulator', 'badread')
        
        if simulator == 'badread' and 'badread' in self.tool_paths:
            # Prepare reference
            ref_to_use = self.config['reference_genome']
            if vcf_file:
                ref_to_use = self.apply_variants_to_reference(sample_id, vcf_file)
                
            cmd = [
                self.tool_paths['badread'],
                'simulate',
                '--reference', ref_to_use,
                '--quantity', f"{coverage}x",
                '--length', str(mean_length), '3000',
                '--error_model', 'nanopore' if platform == 'ont' else 'pacbio',
                '--qscore_model', 'nanopore' if platform == 'ont' else 'pacbio',
            ]
            
            output_file = fastq_dir / f"{sample_id}_long.fastq"
            
            try:
                with open(output_file, 'w') as f:
                    subprocess.run(cmd, stdout=f, check=True)
                    
                # Compress
                subprocess.run(['gzip', str(output_file)], check=True)
                logger.info(f"Long reads generated for {sample_id}")
                
                # Generate BAM
                self.generate_bam(sample_id, 'long')
                
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to generate long reads: {e}")
                
        elif simulator == 'pbsim3' and 'pbsim3' in self.tool_paths:
            self.simulate_with_pbsim(sample_id, vcf_file)
        else:
            logger.error(f"Simulator {simulator} not available")
            
    def apply_variants_to_reference(self, sample_id: str, vcf_file: str) -> str:
        """Apply variants to reference genome"""
        output_dir = Path(self.config['output_directory']) / 'references'
        output_dir.mkdir(exist_ok=True)
        
        output_ref = output_dir / f"{sample_id}_ref.fasta"
        
        # Use bcftools consensus to apply variants
        cmd = [
            'bcftools', 'consensus',
            '-f', self.config['reference_genome'],
            '-o', str(output_ref),
            vcf_file
        ]
        
        try:
            subprocess.run(cmd, check=True)
            return str(output_ref)
        except subprocess.CalledProcessError:
            logger.warning("Failed to apply variants, using original reference")
            return self.config['reference_genome']
            
    def generate_bam(self, sample_id: str, read_type: str):
        """Align reads and generate BAM file"""
        logger.info(f"Generating BAM file for {sample_id}")
        
        output_dir = Path(self.config['output_directory'])
        fastq_dir = output_dir / 'fastq'
        bam_dir = output_dir / 'bam'
        
        ref = self.config['reference_genome']
        
        if read_type == 'short':
            r1 = fastq_dir / f"{sample_id}_R1.fastq.gz"
            r2 = fastq_dir / f"{sample_id}_R2.fastq.gz"
            
            if not (r1.exists() and r2.exists()):
                logger.warning(f"FASTQ files not found for {sample_id}")
                return
                
            # Check for BWA
            if shutil.which('bwa'):
                sam_file = bam_dir / f"{sample_id}.sam"
                
                # BWA alignment
                cmd = [
                    'bwa', 'mem',
                    '-t', str(self.config.get('threads', 4)),
                    ref, str(r1), str(r2)
                ]
                
                try:
                    with open(sam_file, 'w') as f:
                        subprocess.run(cmd, stdout=f, check=True)
                        
                    # Convert to BAM
                    bam_file = bam_dir / f"{sample_id}.bam"
                    subprocess.run([
                        'samtools', 'view',
                        '-b', '-o', str(bam_file),
                        str(sam_file)
                    ], check=True)
                    
                    # Sort and index
                    sorted_bam = bam_dir / f"{sample_id}.sorted.bam"
                    subprocess.run([
                        'samtools', 'sort',
                        '-o', str(sorted_bam),
                        str(bam_file)
                    ], check=True)
                    
                    subprocess.run(['samtools', 'index', str(sorted_bam)], check=True)
                    
                    # Clean up
                    os.remove(sam_file)
                    os.remove(bam_file)
                    
                    logger.info(f"BAM file generated: {sorted_bam}")
                    
                except subprocess.CalledProcessError as e:
                    logger.error(f"Failed to generate BAM: {e}")
            else:
                logger.warning("BWA not found, skipping BAM generation")
                
        elif read_type == 'long':
            reads = fastq_dir / f"{sample_id}_long.fastq.gz"
            
            if not reads.exists():
                logger.warning(f"Long read FASTQ not found for {sample_id}")
                return
                
            # Use minimap2 for long reads
            if shutil.which('minimap2'):
                bam_file = bam_dir / f"{sample_id}_long.bam"
                
                cmd = [
                    'minimap2',
                    '-ax', 'map-ont' if 'ont' in self.config.get('long_read_config', {}).get('platform', 'ont') else 'map-pb',
                    '-t', str(self.config.get('threads', 4)),
                    ref, str(reads)
                ]
                
                try:
                    sam_output = subprocess.run(cmd, capture_output=True, text=True, check=True)
                    
                    # Convert to BAM
                    subprocess.run([
                        'samtools', 'view',
                        '-b', '-o', str(bam_file),
                        '-'
                    ], input=sam_output.stdout, text=True, check=True)
                    
                    # Sort and index
                    sorted_bam = bam_dir / f"{sample_id}_long.sorted.bam"
                    subprocess.run([
                        'samtools', 'sort',
                        '-o', str(sorted_bam),
                        str(bam_file)
                    ], check=True)
                    
                    subprocess.run(['samtools', 'index', str(sorted_bam)], check=True)
                    
                    # Clean up
                    os.remove(bam_file)
                    
                    logger.info(f"Long-read BAM file generated: {sorted_bam}")
                    
                except subprocess.CalledProcessError as e:
                    logger.error(f"Failed to generate long-read BAM: {e}")
            else:
                logger.warning("minimap2 not found, skipping long-read BAM generation")
                
    def generate_readme(self):
        """Generate README file with simulation details"""
        output_dir = Path(self.config['output_directory'])
        readme_path = output_dir / 'README.md'
        
        readme_content = f"""# SeqSimulator Output

## Simulation Date
{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

## Configuration
- Simulation Type: {self.config['simulation_type']}
- Reference Genome: {self.config['reference_genome']}
- Number of Samples: {self.config.get('num_samples', 1)}

## Parameters
"""
        
        if self.config['simulation_type'] in ['short_read', 'hybrid']:
            sr_config = self.config.get('short_read_config', {})
            readme_content += f"""
### Short Read Configuration
- Coverage: {sr_config.get('coverage', 30)}x
- Read Length: {sr_config.get('read_length', 150)}bp
- Insert Size: {sr_config.get('insert_size', 350)}bp
- Simulator: {self.config.get('short_read_simulator', 'art')}
"""
        
        if self.config['simulation_type'] in ['long_read', 'hybrid']:
            lr_config = self.config.get('long_read_config', {})
            readme_content += f"""
### Long Read Configuration
- Coverage: {lr_config.get('coverage', 20)}x
- Mean Length: {lr_config.get('mean_length', 10000)}bp
- Error Rate: {lr_config.get('error_rate', 0.1)}
- Platform: {lr_config.get('platform', 'ont')}
- Simulator: {self.config.get('long_read_simulator', 'badread')}
"""
        
        if 'variant_config' in self.config:
            var_config = self.config['variant_config']
            readme_content += f"""
### Variant Configuration
- SNP Rate: {var_config.get('snp_rate', 0.001)}
- INDEL Rate: {var_config.get('indel_rate', 0.0001)}
- SV Rate: {var_config.get('sv_rate', 0.00001)}
"""
        
        readme_content += f"""
## Output Files

### Directory Structure
```
{self.config['output_directory']}/
├── fastq/       # FASTQ files
├── bam/         # Aligned BAM files
├── vcf/         # Variant files
├── logs/        # Simulation logs
└── reports/     # Quality reports
```

### File Naming Convention
- Short reads: `<sample_id>_R1.fastq.gz`, `<sample_id>_R2.fastq.gz`
- Long reads: `<sample_id>_long.fastq.gz`
- BAM files: `<sample_id>.sorted.bam`
- VCF files: `<sample_id>.vcf.gz`

## Usage with Downstream Tools

### GATK
```bash
gatk HaplotypeCaller \\
    -R {self.config['reference_genome']} \\
    -I bam/<sample>.sorted.bam \\
    -O gatk_output.vcf
```

### Exomiser
```bash
java -jar exomiser.jar \\
    --vcf vcf/<sample>.vcf.gz \\
    --ped <pedigree_file> \\
    --hpo-ids <HPO_terms>
```

### vg (Variation Graph)
```bash
vg construct -r {self.config['reference_genome']} \\
    -v vcf/<sample>.vcf.gz > graph.vg
```

## Contact
For issues or questions, please open an issue on GitHub.
"""
        
        with open(readme_path, 'w') as f:
            f.write(readme_content)
            
        logger.info(f"README generated at {readme_path}")
        
    def run(self):
        """Main execution method"""
        logger.info("Starting simulation...")
        
        num_samples = self.config.get('num_samples', 1)
        sample_prefix = self.config.get('sample_prefix', 'sample')
        
        for i in range(num_samples):
            sample_id = f"{sample_prefix}_{i+1:03d}"
            logger.info(f"Processing sample {sample_id}")
            
            # Generate variants if configured
            vcf_file = None
            if 'variant_config' in self.config:
                vcf_file = self.simulate_variants(sample_id)
                
            # Run simulations based on type
            if self.config['simulation_type'] == 'short_read':
                self.simulate_short_reads(sample_id, vcf_file)
            elif self.config['simulation_type'] == 'long_read':
                self.simulate_long_reads(sample_id, vcf_file)
            elif self.config['simulation_type'] == 'hybrid':
                self.simulate_short_reads(sample_id, vcf_file)
                self.simulate_long_reads(sample_id, vcf_file)
                
        # Generate README
        self.generate_readme()
        
        # Generate summary report
        self.generate_summary_report()
        
        logger.info("Simulation completed successfully!")
        
    def generate_summary_report(self):
        """Generate summary report of simulation"""
        output_dir = Path(self.config['output_directory'])
        report_file = output_dir / 'reports' / 'summary.json'
        
        summary = {
            'timestamp': datetime.now().isoformat(),
            'config': self.config,
            'samples_generated': self.config.get('num_samples', 1),
            'output_files': {}
        }
        
        # List generated files
        for subdir in ['fastq', 'bam', 'vcf']:
            dir_path = output_dir / subdir
            if dir_path.exists():
                files = [f.name for f in dir_path.iterdir() if f.is_file()]
                summary['output_files'][subdir] = files
                
        with open(report_file, 'w') as f:
            json.dump(summary, f, indent=2)
            
        logger.info(f"Summary report generated at {report_file}")


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='SeqSimulator - Simulate short and long-read sequencing data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with configuration file
  python seqsimulator.py -c config.yaml
  
  # Run with configuration and verbose output
  python seqsimulator.py -c config.yaml -v
  
  # Generate example configuration
  python seqsimulator.py --generate-config
        """
    )
    
    parser.add_argument(
        '-c', '--config',
        type=str,
        help='Path to YAML configuration file'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output'
    )
    
    parser.add_argument(
        '--generate-config',
        action='store_true',
        help='Generate example configuration file'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s 1.0.0'
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    if args.generate_config:
        generate_example_config()
        sys.exit(0)
        
    if not args.config:
        parser.error("Configuration file is required. Use -c/--config or --generate-config")
        
    # Run simulator
    simulator = SeqSimulator(args.config)
    simulator.run()


def generate_example_config():
    """Generate example configuration file"""
    example_config = """# SeqSimulator Configuration File
# Generated: {date}

# Simulation type: short_read, long_read, or hybrid
simulation_type: hybrid

# Reference genome (FASTA format, must be indexed)
reference_genome: /path/to/reference/hg38.fa

# Output directory
output_directory: ./simulation_output

# Number of samples to generate
num_samples: 3

# Sample name prefix
sample_prefix: sample

# Number of threads to use
threads: 4

# Chromosomes to simulate (optional, defaults to all)
# Uncomment to simulate specific chromosomes only
# chromosomes:
#   - chr1
#   - chr2
#   - chr22

# Short-read configuration
short_read_config:
  coverage: 30              # Target coverage
  read_length: 150          # Read length in bp
  insert_size: 350          # Insert size for paired-end
  insert_sd: 50             # Insert size standard deviation
  
# Short-read simulator: art, dwgsim
short_read_simulator: art

# Long-read configuration
long_read_config:
  coverage: 20              # Target coverage
  mean_length: 10000        # Mean read length
  error_rate: 0.1           # Base error rate
  platform: ont             # Platform: ont or pacbio
  
# Long-read simulator: badread, pbsim
long_read_simulator: badread

# Variant configuration
variant_config:
  # Option 1: Use existing population VCF
  # population_vcf: /path/to/population.vcf.gz
  
  # Option 2: Generate random variants
  snp_rate: 0.001           # SNPs per base
  indel_rate: 0.0001        # INDELs per base
  sv_rate: 0.00001          # Structural variants per base
  
  # Variant allele frequency distribution
  af_distribution:
    rare: 0.7               # Proportion of rare variants (AF < 0.01)
    low: 0.2                # Proportion of low frequency (0.01 < AF < 0.05)
    common: 0.1             # Proportion of common variants (AF > 0.05)

# Advanced options
advanced:
  keep_intermediate: false   # Keep intermediate files
  seed: 42                  # Random seed for reproducibility
  quality_model: default    # Quality score model
""".format(date=datetime.now().strftime('%Y-%m-%d'))
    
    config_file = 'example_config.yaml'
    with open(config_file, 'w') as f:
        f.write(example_config)
        
    print(f"Example configuration file generated: {config_file}")
    print("Please edit this file with your specific parameters and paths.")


if __name__ == '__main__':
    main()
