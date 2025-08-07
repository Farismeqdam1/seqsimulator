#!/usr/bin/env python3
"""
Validate SeqSimulator output files
"""

import os
import sys
import gzip
import json
import argparse
from pathlib import Path
import pysam
from Bio import SeqIO
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


class OutputValidator:
    """Validate simulation output files"""
    
    def __init__(self, output_dir):
        self.output_dir = Path(output_dir)
        self.errors = []
        self.warnings = []
        self.stats = {}
        
    def validate_all(self):
        """Run all validation checks"""
        logger.info(f"Validating output directory: {self.output_dir}")
        
        # Check directory structure
        self.check_directory_structure()
        
        # Validate FASTQ files
        self.validate_fastq_files()
        
        # Validate BAM files
        self.validate_bam_files()
        
        # Validate VCF files
        self.validate_vcf_files()
        
        # Check README and reports
        self.check_documentation()
        
        # Generate validation report
        self.generate_report()
        
        return len(self.errors) == 0
        
    def check_directory_structure(self):
        """Check if all required directories exist"""
        required_dirs = ['fastq', 'bam', 'vcf', 'logs', 'reports']
        
        for dir_name in required_dirs:
            dir_path = self.output_dir / dir_name
            if not dir_path.exists():
                self.errors.append(f"Missing directory: {dir_name}")
            else:
                logger.info(f"✓ Directory exists: {dir_name}")
                
    def validate_fastq_files(self):
        """Validate FASTQ files"""
        fastq_dir = self.output_dir / 'fastq'
        
        if not fastq_dir.exists():
            return
            
        fastq_files = list(fastq_dir.glob('*.fastq.gz')) + list(fastq_dir.glob('*.fq.gz'))
        
        if not fastq_files:
            self.warnings.append("No FASTQ files found")
            return
            
        logger.info(f"Found {len(fastq_files)} FASTQ files")
        
        for fastq_file in fastq_files:
            try:
                # Check if file is not empty
                if fastq_file.stat().st_size == 0:
                    self.errors.append(f"Empty FASTQ file: {fastq_file.name}")
                    continue
                    
                # Parse and validate FASTQ format
                read_count = 0
                total_length = 0
                min_qual = float('inf')
                max_qual = 0
                
                with gzip.open(fastq_file, 'rt') as handle:
                    for record in SeqIO.parse(handle, 'fastq'):
                        read_count += 1
                        total_length += len(record.seq)
                        
                        # Check quality scores
                        quals = record.letter_annotations.get('phred_quality', [])
                        if quals:
                            min_qual = min(min_qual, min(quals))
                            max_qual = max(max_qual, max(quals))
                            
                        # Validate first 100 reads only for speed
                        if read_count >= 100:
                            break
                            
                if read_count > 0:
                    avg_length = total_length / read_count
                    self.stats[fastq_file.name] = {
                        'reads_checked': read_count,
                        'avg_length': avg_length,
                        'min_qual': min_qual if min_qual != float('inf') else 0,
                        'max_qual': max_qual
                    }
                    logger.info(f"✓ {fastq_file.name}: {read_count} reads validated, avg length: {avg_length:.1f}")
                else:
                    self.errors.append(f"No valid reads in {fastq_file.name}")
                    
            except Exception as e:
                self.errors.append(f"Error validating {fastq_file.name}: {str(e)}")
                
    def validate_bam_files(self):
        """Validate BAM files"""
        bam_dir = self.output_dir / 'bam'
        
        if not bam_dir.exists():
            return
            
        bam_files = list(bam_dir.glob('*.bam'))
        
        if not bam_files:
            self.warnings.append("No BAM files found")
            return
            
        logger.info(f"Found {len(bam_files)} BAM files")
        
        for bam_file in bam_files:
            try:
                # Open BAM file
                bamfile = pysam.AlignmentFile(str(bam_file), 'rb')
                
                # Check if index exists
                bai_file = Path(str(bam_file) + '.bai')
                if not bai_file.exists():
                    self.warnings.append(f"Missing BAM index: {bai_file.name}")
                    
                # Get statistics
                total_reads = 0
                mapped_reads = 0
                
                for read in bamfile.head(1000):  # Check first 1000 reads
                    total_reads += 1
                    if not read.is_unmapped:
                        mapped_reads += 1
                        
                if total_reads > 0:
                    mapping_rate = (mapped_reads / total_reads) * 100
                    self.stats[bam_file.name] = {
                        'reads_checked': total_reads,
                        'mapped_reads': mapped_reads,
                        'mapping_rate': mapping_rate
                    }
                    logger.info(f"✓ {bam_file.name}: {mapping_rate:.1f}% mapping rate")
                    
                bamfile.close()
                
            except Exception as e:
                self.errors.append(f"Error validating {bam_file.name}: {str(e)}")
                
    def validate_vcf_files(self):
        """Validate VCF files"""
        vcf_dir = self.output_dir / 'vcf'
        
        if not vcf_dir.exists():
            return
            
        vcf_files = list(vcf_dir.glob('*.vcf.gz')) + list(vcf_dir.glob('*.vcf'))
        
        if not vcf_files:
            self.warnings.append("No VCF files found")
            return
            
        logger.info(f"Found {len(vcf_files)} VCF files")
        
        for vcf_file in vcf_files:
            try:
                # Open VCF file
                vcf = pysam.VariantFile(str(vcf_file))
                
                # Check header
                if not vcf.header:
                    self.errors.append(f"Missing VCF header: {vcf_file.name}")
                    continue
                    
                # Count variants by type
                snp_count = 0
                indel_count = 0
                sv_count = 0
                total_variants = 0
                
                for variant in vcf.fetch():
                    total_variants += 1
                    
                    # Classify variant type
                    ref_len = len(variant.ref)
                    alt_lens = [len(alt) for alt in variant.alts] if variant.alts else [0]
                    
                    if ref_len == 1 and all(alt_len == 1 for alt_len in alt_lens):
                        snp_count += 1
                    elif ref_len < 50 or any(alt_len < 50 for alt_len in alt_lens):
                        indel_count += 1
                    else:
                        sv_count += 1
                        
                    # Check first 1000 variants only
                    if total_variants >= 1000:
                        break
                        
                self.stats[vcf_file.name] = {
                    'total_variants': total_variants,
                    'snps': snp_count,
                    'indels': indel_count,
                    'svs': sv_count
                }
                
                logger.info(f"✓ {vcf_file.name}: {total_variants} variants (SNPs: {snp_count}, INDELs: {indel_count}, SVs: {sv_count})")
                
                vcf.close()
                
            except Exception as e:
                self.errors.append(f"Error validating {vcf_file.name}: {str(e)}")
                
    def check_documentation(self):
        """Check for README and report files"""
        readme_file = self.output_dir / 'README.md'
        if not readme_file.exists():
            self.warnings.append("README.md not found")
        else:
            logger.info("✓ README.md found")
            
        summary_file = self.output_dir / 'reports' / 'summary.json'
        if summary_file.exists():
            try:
                with open(summary_file, 'r') as f:
                    summary = json.load(f)
                logger.info("✓ Summary report found and valid")
            except Exception as e:
                self.errors.append(f"Invalid summary.json: {str(e)}")
        else:
            self.warnings.append("summary.json not found")
            
    def generate_report(self):
        """Generate validation report"""
        report = {
            'validation_passed': len(self.errors) == 0,
            'errors': self.errors,
            'warnings': self.warnings,
            'statistics': self.stats
        }
        
        report_file = self.output_dir / 'validation_report.json'
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
            
        logger.info(f"\nValidation Report saved to: {report_file}")
        
        # Print summary
        print("\n" + "="*50)
        print("VALIDATION SUMMARY")
        print("="*50)
        
        if len(self.errors) == 0:
            print("✓ Validation PASSED")
        else:
            print("✗ Validation FAILED")
            print(f"\nErrors ({len(self.errors)}):")
            for error in self.errors:
                print(f"  - {error}")
                
        if self.warnings:
            print(f"\nWarnings ({len(self.warnings)}):")
            for warning in self.warnings:
                print(f"  - {warning}")
                
        print("\n" + "="*50)


def main():
    parser = argparse.ArgumentParser(description='Validate SeqSimulator output')
    parser.add_argument(
        'output_dir',
        help='Path to SeqSimulator output directory'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Verbose output'
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    if not os.path.exists(args.output_dir):
        logger.error(f"Output directory not found: {args.output_dir}")
        sys.exit(1)
        
    validator = OutputValidator(args.output_dir)
    success = validator.validate_all()
    
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
