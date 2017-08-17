/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package de.helmholtz_muenchen.ibis.ngs.pindel;

import javax.swing.JFileChooser;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentDate;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentLabel;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentNumberEdit;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDate;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.BAMCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;

/**
 * <code>NodeDialog</code> for the "Pindel" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more
 * complex dialog please derive directly from
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Marie-Sophie Friedl
 * @author Maximilian Hastreiter
 */
public class PindelNodeDialog extends HTExecutorNodeDialog {
	
	public PindelNodeDialog(){
		super(BAMCell.TYPE.getPreferredValueClass(), FileCell.TYPE.getPreferredValueClass());
	}

	public void addToolDialogComponents() {

		final SettingsModelString pindel = new SettingsModelString(PindelNodeModel.CFGKEY_PINDEL,
				PindelNodeModel.DEF_PINDEL);
		;
		final SettingsModelString refseqfile = new SettingsModelString(PindelNodeModel.CFGKEY_REFSEQFILE, "");
		final SettingsModelBoolean interval = new SettingsModelBoolean(PindelNodeModel.CFGKEY_INTERVAL,
				PindelNodeModel.DEF_INTERVAL);
		;
		final SettingsModelString chrom = new SettingsModelString(PindelNodeModel.CFGKEY_CHROM,
				PindelNodeModel.DEF_CHROM);
		final SettingsModelIntegerBounded start = new SettingsModelIntegerBounded(PindelNodeModel.CFGKEY_START,
				PindelNodeModel.DEF_START, PindelNodeModel.MIN_START, PindelNodeModel.MAX_START);
		final SettingsModelIntegerBounded end = new SettingsModelIntegerBounded(PindelNodeModel.CFGKEY_END,
				PindelNodeModel.DEF_END, PindelNodeModel.MIN_END, PindelNodeModel.MAX_END);
		final SettingsModelString config_file = new SettingsModelString(PindelNodeModel.CFGKEY_CONFIG_FILE,
				PindelNodeModel.DEF_CONFIG_FILE);
		final SettingsModelBoolean create_config = new SettingsModelBoolean(PindelNodeModel.CFGKEY_CREATE_CONFIG,
				PindelNodeModel.DEF_CREATE_CONFIG);
		final SettingsModelBoolean vcf_out = new SettingsModelBoolean(PindelNodeModel.CFGKEY_VCF_OUT,
				PindelNodeModel.DEF_VCF_OUT);
		final SettingsModelString pindel2vcf = new SettingsModelString(PindelNodeModel.CFGKEY_VCF2PINDEL,
				PindelNodeModel.DEF_VCF2PINDEL);
		final SettingsModelIntegerBounded threads = new SettingsModelIntegerBounded(PindelNodeModel.CFGKEY_THREADS,
				PindelNodeModel.DEF_THREADS, PindelNodeModel.MIN_THREADS, PindelNodeModel.MAX_THREADS);
		final SettingsModelIntegerBounded bin_size = new SettingsModelIntegerBounded(PindelNodeModel.CFGKEY_BIN_SIZE,
				PindelNodeModel.DEF_BIN_SIZE, PindelNodeModel.MIN_BIN_SIZE, PindelNodeModel.MAX_BIN_SIZE);

		// pindel params
		final SettingsModelInteger min_match_bases = new SettingsModelIntegerBounded(
				PindelNodeModel.CFGKEY_MIN_MATCH_BASES, PindelNodeModel.DEF_MIN_MATCH_BASES,
				PindelNodeModel.MIN_MIN_MATCH_BASES, PindelNodeModel.MAX_MIN_MATCH_BASES);
		final SettingsModelInteger additional_mismatch = new SettingsModelIntegerBounded(
				PindelNodeModel.CFGKEY_ADDITIONAL_MISMATCH, PindelNodeModel.DEF_ADDITIONAL_MISMATCH,
				PindelNodeModel.MIN_ADDITIONAL_MISMATCH, PindelNodeModel.MAX_ADDITIONAL_MISMATCH);
		final SettingsModelInteger min_match_breakpoint = new SettingsModelIntegerBounded(
				PindelNodeModel.CFGKEY_MIN_MATCH_BP, PindelNodeModel.DEF_MIN_MATCH_BP, PindelNodeModel.MIN_MIN_MATCH_BP,
				PindelNodeModel.MAX_MIN_MATCH_BP);
		final SettingsModelDoubleBounded seq_err = new SettingsModelDoubleBounded(PindelNodeModel.CFGKEY_SEQ_ERR,
				PindelNodeModel.DEF_SEQ_ERROR, PindelNodeModel.MIN_SEQ_ERROR, PindelNodeModel.MAX_SEQ_ERROR);
		final SettingsModelDoubleBounded max_mismatch_rate = new SettingsModelDoubleBounded(
				PindelNodeModel.CFGKEY_MAX_MISMATCH_RATE, PindelNodeModel.DEF_MAX_MISMATCH_RATE,
				PindelNodeModel.MIN_MAX_MISMATCH_RATE, PindelNodeModel.MAX_MAX_MISMATCH_RATE);

		// pindel2vcf params
		final SettingsModelBoolean use_ref_filename = new SettingsModelBoolean(PindelNodeModel.CFGKEY_USE_REF_FILENAME,
				PindelNodeModel.DEF_USE_REF_FILENAME);
		final SettingsModelString refname = new SettingsModelString(PindelNodeModel.CFGKEY_REFNAME,
				PindelNodeModel.DEF_REFNAME);
		final SettingsModelBoolean use_cur_date = new SettingsModelBoolean(PindelNodeModel.CFGKEY_USE_CUR_DATE,
				PindelNodeModel.DEF_USE_CUR_DATE);
		final SettingsModelDate refdate = new SettingsModelDate(PindelNodeModel.CFGKEY_REFDATE);
		final SettingsModelIntegerBounded min_reads = new SettingsModelIntegerBounded(PindelNodeModel.CFGKEY_MIN_READS,
				PindelNodeModel.DEF_MIN_READS, PindelNodeModel.MIN_MIN_READS, PindelNodeModel.MAX_MIN_READS);
		final SettingsModelDoubleBounded hetero_frac = new SettingsModelDoubleBounded(
				PindelNodeModel.CFGKEY_HETERO_FRAC, PindelNodeModel.DEF_HETERO_FRAC, PindelNodeModel.MIN_HETERO_FRAC,
				PindelNodeModel.MAX_HETERO_FRAC);
		final SettingsModelDoubleBounded homo_frac = new SettingsModelDoubleBounded(PindelNodeModel.CFGKEY_HOMO_FRAC,
				PindelNodeModel.DEF_HOMO_FRAC, PindelNodeModel.MIN_HOMO_FRAC, PindelNodeModel.MAX_HOMO_FRAC);
		final SettingsModelBoolean gatk_comp = new SettingsModelBoolean(PindelNodeModel.CFGKEY_GATK_COMP,
				PindelNodeModel.DEF_GATK_COMP);

		final SettingsModelIntegerBounded min_supp_reads = new SettingsModelIntegerBounded(
				PindelNodeModel.CFGKEY_MIN_SUPP_READS, PindelNodeModel.DEF_MIN_SUPP_READS,
				PindelNodeModel.MIN_MIN_SUPP_READS, PindelNodeModel.MAX_MIN_SUPP_READS);
		final SettingsModelBoolean both_strands = new SettingsModelBoolean(PindelNodeModel.CFGKEY_BOTH_STRANDS,
				PindelNodeModel.DEF_BOTH_STRANDS);
		final SettingsModelIntegerBounded min_size = new SettingsModelIntegerBounded(PindelNodeModel.CFGKEY_MIN_SIZE,
				PindelNodeModel.DEF_MIN_SIZE, PindelNodeModel.MIN_MIN_SIZE, PindelNodeModel.MAX_MIN_SIZE);
		final SettingsModelBoolean limit_size = new SettingsModelBoolean(PindelNodeModel.CFGKEY_LIMIT_SIZE,
				PindelNodeModel.DEF_LIMIT_SIZE);
		final SettingsModelIntegerBounded max_size = new SettingsModelIntegerBounded(PindelNodeModel.CFGKEY_MAX_SIZE,
				PindelNodeModel.DEF_MAX_SIZE, PindelNodeModel.MIN_MAX_SIZE, PindelNodeModel.MAX_MAX_SIZE);

		addPrefPageSetting(pindel, IBISKNIMENodesPlugin.PINDEL);
		addPrefPageSetting(pindel2vcf, IBISKNIMENodesPlugin.PINDEL2VCF);
		addPrefPageSetting(refseqfile, IBISKNIMENodesPlugin.REF_GENOME);

		createNewGroup("Interval for variant calling");
		addDialogComponent(
				new DialogComponentBoolean(interval, "Restrict variant calling to a certain genomic region"));
		setHorizontalPlacement(true);
		addDialogComponent(new DialogComponentString(chrom, "Chromosome"));
		addDialogComponent(new DialogComponentNumberEdit(start, "Start", 8));
		addDialogComponent(new DialogComponentNumberEdit(end, "End", 8));
		chrom.setEnabled(false);
		start.setEnabled(false);
		end.setEnabled(false);
		setHorizontalPlacement(false);
		addDialogComponent(new DialogComponentLabel("(Chromosome name has to match reference and BAM file header)"));

		// text/number fields for chromosome, start and end only active if
		// interval is used
		interval.addChangeListener(new ChangeListener() {

			@Override
			public void stateChanged(ChangeEvent e) {
				chrom.setEnabled(interval.getBooleanValue());
				start.setEnabled(interval.getBooleanValue());
				end.setEnabled(interval.getBooleanValue());
			}
		});

		createNewGroup("Path to Pindel config file");
		addDialogComponent(new DialogComponentFileChooser(config_file, "pindelconfig", JFileChooser.OPEN_DIALOG, false,
				".config"));
		create_config.setEnabled(true);
		addDialogComponent(new DialogComponentBoolean(create_config,
				"Create config file (requires PicardTools: CollectInsertMetrics as previous node)"));

		// disables file chooser if config file is create by this node
		create_config.addChangeListener(new ChangeListener() {

			@Override
			public void stateChanged(ChangeEvent e) {
				config_file.setEnabled(!create_config.getBooleanValue());
			}
		});

		createNewGroup("Output");
		addDialogComponent(new DialogComponentBoolean(vcf_out,
				"Convert Pindel output (deletions and small insertions) to VCF format"));

		createNewGroup("Runtime and memeroy usage");
		setHorizontalPlacement(true);
		addDialogComponent(new DialogComponentNumber(threads, "Number of threads", 1, 5));
		addDialogComponent(new DialogComponentNumber(bin_size, "Bin size", 1, 5));
		addDialogComponent(new DialogComponentLabel("(Mb of reference in memory)"));
		setHorizontalPlacement(false);

		// }

		// private void PindelParams(){

		createNewTab("Pindel Parameters");

		createNewGroup("Minimum number of matching bases");
		addDialogComponent(new DialogComponentLabel(
				"Only consider reads as evidence if they map with more than this number of bases:"));
		addDialogComponent(new DialogComponentNumber(min_match_bases, "Mismatching bases per read", 1, 5));

		createNewGroup("Mismatch threshold");
		addDialogComponent(new DialogComponentLabel(
				"Do not align a read if there is another mapping position below this threshold:"));
		addDialogComponent(new DialogComponentNumber(additional_mismatch, "Mismatching bases per alignment", 1, 5));

		createNewGroup("Number of perfect matches at breakpoints");
		addDialogComponent(
				new DialogComponentLabel("Number of perfectly matching bases around a breakpoint of a split read:"));
		addDialogComponent(new DialogComponentNumber(min_match_breakpoint, "Number of perfect matches", 1, 5));

		createNewGroup("Sequencing error rate");
		addDialogComponent(new DialogComponentLabel("Expected fraction of sequencing errors:"));
		addDialogComponent(new DialogComponentNumber(seq_err, "Error rate", 0.01, 5));
		createNewGroup("Maximum allowed mismatch rate");

		addDialogComponent(
				new DialogComponentLabel("Only consider aligned reads with mismatch rate below this fraction:"));
		addDialogComponent(
				new DialogComponentNumber(max_mismatch_rate, "Fraction of mismatching bases per read", 0.01, 5));

		// }

		// private void Pindel2VCFParams(){

		createNewTab("Pindel2vcf Parameters");

		createNewGroup("Reference sequence");
		addDialogComponent(new DialogComponentBoolean(use_ref_filename, "Use file name as reference name"));
		addDialogComponent(new DialogComponentString(refname, "Name of the reference sequence:", true, 10));
		refname.setEnabled(false);

		use_ref_filename.addChangeListener(new ChangeListener() {

			@Override
			public void stateChanged(ChangeEvent e) {

				refname.setEnabled(!use_ref_filename.getBooleanValue());

			}
		});

		addDialogComponent(new DialogComponentBoolean(use_cur_date, "Use current date"));
		addDialogComponent(new DialogComponentDate(refdate, "Date of the version of the reference sequence"));
		refdate.setEnabled(false);

		use_cur_date.addChangeListener(new ChangeListener() {

			@Override
			public void stateChanged(ChangeEvent e) {

				refdate.setEnabled(!use_cur_date.getBooleanValue());

			}
		});

		createNewGroup("Genotype");
		addDialogComponent(new DialogComponentNumber(min_reads, "Minimum reads to report genotype", 1, 5));
		addDialogComponent(
				new DialogComponentNumber(hetero_frac, "Proportion of reads defined as heterozygous", 0.01, 5));
		addDialogComponent(new DialogComponentNumber(homo_frac, "Proportion of reads defined as homozygous", 0.01, 5));
		addDialogComponent(new DialogComponentBoolean(gatk_comp, "Output GATK-compatible genotypes (recommended)"));

		createNewGroup("Filter");
		addDialogComponent(new DialogComponentBoolean(both_strands,
				"Only output variants that are supported by reads on both strands"));
		addDialogComponent(new DialogComponentNumber(min_supp_reads, "Minimum number of supporting reads", 1, 5));
		addDialogComponent(new DialogComponentNumber(min_size, "Minimum variant size", 1, 5));
		setHorizontalPlacement(true);
		addDialogComponent(new DialogComponentBoolean(limit_size, "Limit variant size"));
		addDialogComponent(new DialogComponentNumber(max_size, "Maximum variant size", 1, 5));
		max_size.setEnabled(false);
		setHorizontalPlacement(false);

		limit_size.addChangeListener(new ChangeListener() {

			@Override
			public void stateChanged(ChangeEvent e) {

				max_size.setEnabled(limit_size.getBooleanValue());
			}
		});
	}

}
