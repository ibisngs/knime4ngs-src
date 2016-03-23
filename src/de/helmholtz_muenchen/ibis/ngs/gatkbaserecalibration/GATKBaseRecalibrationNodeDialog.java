package de.helmholtz_muenchen.ibis.ngs.gatkbaserecalibration;

import javax.swing.JFileChooser;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

public class GATKBaseRecalibrationNodeDialog extends HTExecutorNodeDialog {

	public void addToolDialogComponents() {

		final SettingsModelString gatk = new SettingsModelString(GATKBaseRecalibrationNodeModel.CFGKEY_GATK, "");
		final SettingsModelString reffile = new SettingsModelString(GATKBaseRecalibrationNodeModel.CFGKEY_REF_GENOME,
				"");
		final SettingsModelBoolean use_phase1_1000G = new SettingsModelBoolean(
				GATKBaseRecalibrationNodeModel.CFGKEY_USE_PHASE1_1000G,
				GATKBaseRecalibrationNodeModel.DEF_USE_PHASE1_1000G);
		final SettingsModelString phase1_1000G_file = new SettingsModelString(
				GATKBaseRecalibrationNodeModel.CFGKEY_PHASE1_1000G_FILE,
				GATKBaseRecalibrationNodeModel.DEF_PHASE1_1000G_FILE);
		final SettingsModelBoolean use_mills_1000G = new SettingsModelBoolean(
				GATKBaseRecalibrationNodeModel.CFGKEY_USE_MILLS_1000G,
				GATKBaseRecalibrationNodeModel.DEF_USE_MILLS_1000G);
		final SettingsModelString mills_1000G_file = new SettingsModelString(
				GATKBaseRecalibrationNodeModel.CFGKEY_MILLS_1000G_FILE,
				GATKBaseRecalibrationNodeModel.DEF_MILLS_1000G_FILE);
		final SettingsModelBoolean use_dbsnp = new SettingsModelBoolean(GATKBaseRecalibrationNodeModel.CFGKEY_USE_DBSNP,
				GATKBaseRecalibrationNodeModel.DEF_USE_DBSNP);
		final SettingsModelString dbsnp_file = new SettingsModelString(GATKBaseRecalibrationNodeModel.CFGKEY_DBSNP_FILE,
				GATKBaseRecalibrationNodeModel.DEF_DBSNP_FILE);
		final SettingsModelBoolean use_interval = new SettingsModelBoolean(
				GATKBaseRecalibrationNodeModel.CFGKEY_USE_INTERVAL, GATKBaseRecalibrationNodeModel.DEF_USE_INTERVAL);
		final SettingsModelString interval_file = new SettingsModelString(
				GATKBaseRecalibrationNodeModel.CFGKEY_INTERVAL_FILE, GATKBaseRecalibrationNodeModel.DEF_INTERVAL_FILE);
		final SettingsModelBoolean create_plots = new SettingsModelBoolean(
				GATKBaseRecalibrationNodeModel.CFGKEY_CREATE_PLOTS, GATKBaseRecalibrationNodeModel.DEF_CREATE_PLOTS);
		final SettingsModelIntegerBounded cpu_threads = new SettingsModelIntegerBounded(
				GATKBaseRecalibrationNodeModel.CFGKEY_CPU_THREADS, GATKBaseRecalibrationNodeModel.DEF_CPU_THREADS,
				GATKBaseRecalibrationNodeModel.MIN_CPU_THREADS, GATKBaseRecalibrationNodeModel.MAX_CPU_THREADS);
		final SettingsModelIntegerBounded memory_usage = new SettingsModelIntegerBounded(
				GATKBaseRecalibrationNodeModel.CFGKEY_JAVAMEMORY, GATKBaseRecalibrationNodeModel.DEF_NUM_JAVAMEMORY,
				GATKBaseRecalibrationNodeModel.MIN_NUM_JAVAMEMORY, GATKBaseRecalibrationNodeModel.MAX_NUM_JAVAMEMORY);

//		final SettingsModelBoolean context_cov = new SettingsModelBoolean(
//				GATKBaseRecalibrationNodeModel.CFGKEY_CONTEXT_COV, GATKBaseRecalibrationNodeModel.DEF_CONTEXT_COV);
//		final SettingsModelBoolean cycle_cov = new SettingsModelBoolean(GATKBaseRecalibrationNodeModel.CFGKEY_CYCLE_COV,
//				GATKBaseRecalibrationNodeModel.DEF_CYCLE_COV);
//		final SettingsModelBoolean rep_len_cov = new SettingsModelBoolean(
//				GATKBaseRecalibrationNodeModel.CFGKEY_REP_LEN_COV, GATKBaseRecalibrationNodeModel.DEF_REP_LEN_COV);
//		final SettingsModelBoolean rep_unit_cov = new SettingsModelBoolean(
//				GATKBaseRecalibrationNodeModel.CFGKEY_REP_UNIT_COV, GATKBaseRecalibrationNodeModel.DEF_REP_UNIT_COV);
		final SettingsModelIntegerBounded low_qual_tail = new SettingsModelIntegerBounded(
				GATKBaseRecalibrationNodeModel.CFGKEY_LOW_QUAL_TAIL, GATKBaseRecalibrationNodeModel.DEF_LOW_QUAL_TAIL,
				GATKBaseRecalibrationNodeModel.MIN_LOW_QUAL_TAIL, GATKBaseRecalibrationNodeModel.MAX_LOW_QUAL_TAIL);
		final SettingsModelDoubleBounded gap_open = new SettingsModelDoubleBounded(
				GATKBaseRecalibrationNodeModel.CFGKEY_GAP_OPEN, GATKBaseRecalibrationNodeModel.DEF_GAP_OPEN,
				GATKBaseRecalibrationNodeModel.MIN_GAP_OPEN, GATKBaseRecalibrationNodeModel.MAX_GAP_OPEN);
		final SettingsModelIntegerBounded deletion_def_qual = new SettingsModelIntegerBounded(
				GATKBaseRecalibrationNodeModel.CFGKEY_DELETION_DEF_QUAL,
				GATKBaseRecalibrationNodeModel.DEF_DELETION_DEF_QUAL,
				GATKBaseRecalibrationNodeModel.MIN_DELETION_DEF_QUAL,
				GATKBaseRecalibrationNodeModel.MAX_DELETION_DEF_QUAL);
		final SettingsModelIntegerBounded insertion_def_qual = new SettingsModelIntegerBounded(
				GATKBaseRecalibrationNodeModel.CFGKEY_INSERTION_DEF_QUAL,
				GATKBaseRecalibrationNodeModel.DEF_INSERTION_DEF_QUAL,
				GATKBaseRecalibrationNodeModel.MIN_INSERTION_DEF_QUAL,
				GATKBaseRecalibrationNodeModel.MAX_INSERTION_DEF_QUAL);
		final SettingsModelIntegerBounded indel_context_size = new SettingsModelIntegerBounded(
				GATKBaseRecalibrationNodeModel.CFGKEY_INDEL_CONTEXT_SIZE,
				GATKBaseRecalibrationNodeModel.DEF_INDEL_CONTEXT_SIZE,
				GATKBaseRecalibrationNodeModel.MIN_INDEL_CONTEXT_SIZE,
				GATKBaseRecalibrationNodeModel.MAX_INDEL_CONTEXT_SIZE);
		final SettingsModelIntegerBounded mismatch_def_qual = new SettingsModelIntegerBounded(
				GATKBaseRecalibrationNodeModel.CFGKEY_MISMATCH_DEF_QUAL,
				GATKBaseRecalibrationNodeModel.DEF_MISMATCH_DEF_QUAL,
				GATKBaseRecalibrationNodeModel.MIN_MISMATCH_DEF_QUAL,
				GATKBaseRecalibrationNodeModel.MAX_MISMACTH_DEF_QUAL);
		final SettingsModelIntegerBounded mismatch_context_size = new SettingsModelIntegerBounded(
				GATKBaseRecalibrationNodeModel.CFGKEY_MISMATCH_CONTEXT_SIZE,
				GATKBaseRecalibrationNodeModel.DEF_MISMATCH_CONTEXT_SIZE,
				GATKBaseRecalibrationNodeModel.MIN_MISMATCH_CONTEXT_SIZE,
				GATKBaseRecalibrationNodeModel.MAX_MISMATCH_CONTEXT_SIZE);
		final SettingsModelIntegerBounded max_cycles = new SettingsModelIntegerBounded(
				GATKBaseRecalibrationNodeModel.CFGKEY_MAX_CYCLES, GATKBaseRecalibrationNodeModel.DEF_MAX_CYCLES,
				GATKBaseRecalibrationNodeModel.MIN_MAX_CYCLES, GATKBaseRecalibrationNodeModel.MAX_MAX_CYCLES);
		final SettingsModelBoolean simplify_out = new SettingsModelBoolean(
				GATKBaseRecalibrationNodeModel.CFGKEY_SIMPLIFY_OUT, GATKBaseRecalibrationNodeModel.DEF_SIMPLIY_OUT);
		final SettingsModelOptionalString m_br_opt_flags = new SettingsModelOptionalString(
				GATKBaseRecalibrationNodeModel.CFGKEY_BR_OPT_FLAGS, "", false);
		final SettingsModelOptionalString m_ac_opt_flags = new SettingsModelOptionalString(
				GATKBaseRecalibrationNodeModel.CFGKEY_AC_OPT_FLAGS, "", false);
		final SettingsModelOptionalString m_pr_opt_flags = new SettingsModelOptionalString(
				GATKBaseRecalibrationNodeModel.CFGKEY_PR_OPT_FLAGS, "", false);

		// Proxy options
		// private final SettingsModelBoolean useproxy = new
		// SettingsModelBoolean(GATKBaseRecalibrationNodeModel.CFGKEY_USEPROXY,
		// false);
		// final SettingsModelString proxyhost = new
		// SettingsModelString(GATKBaseRecalibrationNodeModel.CFGKEY_PROXYHOST,
		// null);
		// final SettingsModelString proxyport = new
		// SettingsModelString(GATKBaseRecalibrationNodeModel.CFGKEY_PROXYPORT,
		// null);
		// private final SettingsModelBoolean useproxyauth = new
		// SettingsModelBoolean(GATKBaseRecalibrationNodeModel.CFGKEY_USEPROXYAUTH,
		// false);
		// final SettingsModelString proxyuser = new
		// SettingsModelString(GATKBaseRecalibrationNodeModel.CFGKEY_PROXYUSER,
		// null);
		// final SettingsModelString proxypassword = new
		// SettingsModelString(GATKBaseRecalibrationNodeModel.CFGKEY_PROXYPASSWORD,
		// null);

		addPrefPageSetting(gatk, IBISKNIMENodesPlugin.GATK);
		addPrefPageSetting(reffile, IBISKNIMENodesPlugin.REF_GENOME);
		addPrefPageSetting(phase1_1000G_file, IBISKNIMENodesPlugin.RES_1000G_INDELS);
		addPrefPageSetting(mills_1000G_file, IBISKNIMENodesPlugin.RES_MILLS);
		addPrefPageSetting(dbsnp_file, IBISKNIMENodesPlugin.RES_DBSNP);

		createNewGroup("Sets of known polymorphisms (at least one set has to be chosen)");
		addDialogComponent(new DialogComponentBoolean(use_phase1_1000G, "Use 1000 genomes phase1 indel set"));
		addDialogComponent(
				new DialogComponentBoolean(use_mills_1000G, "Use Mills and 1000 genomes gold standard indel set"));
		addDialogComponent(new DialogComponentBoolean(use_dbsnp, "Use dbSNP variant set"));

		// interval file for realignment
		createNewGroup("Interval for recalibration");
		addDialogComponent(
				new DialogComponentBoolean(use_interval, "Restrict recalibration to certain genomic regions"));
		interval_file.setEnabled(false);
		addDialogComponent(new DialogComponentFileChooser(interval_file, "ifile2", JFileChooser.OPEN_DIALOG, false,
				".bed", ".intervals"));

		// add change listener that enables file chooser for interval file
		use_interval.addChangeListener(new ChangeListener() {

			public void stateChanged(ChangeEvent e) {

				interval_file.setEnabled(use_interval.getBooleanValue());
			}

		});

		createNewGroup("AnalyzeCovariates");
		addDialogComponent(new DialogComponentBoolean(create_plots, "Create before/after plots for base qualities"));
		addDialogComponent(new DialogComponentOptionalString(m_ac_opt_flags,"Optional flags"));
		
		createNewGroup("PrintReads");
		addDialogComponent(new DialogComponentBoolean(simplify_out,"Remove all additional tags from BAM file (except read group)"));
		addDialogComponent(new DialogComponentOptionalString(m_pr_opt_flags,"Optional flags"));
		
		createNewGroup("General Options");
		addDialogComponent(new DialogComponentNumber(cpu_threads, "CPU threads", 1, 5));
		addDialogComponent(new DialogComponentNumber(memory_usage, "Shared Java Memory (GB) for all CPU threads", 1, 5));

		createNewTab("BaseRecalibrator");
//		createNewGroup("Calculation of covariates");
//		addDialogComponent(new DialogComponentLabel("Quality score covariate"));
//		addDialogComponent(new DialogComponentLabel("Read group covariate"));
//		addDialogComponent(new DialogComponentBoolean(context_cov, "Context covariate"));
//		addDialogComponent(new DialogComponentBoolean(cycle_cov, "Cycle covariate"));
//		addDialogComponent(new DialogComponentBoolean(rep_len_cov, "Repeat length covariate"));
//		addDialogComponent(new DialogComponentBoolean(rep_unit_cov, "Repeat unit covariate"));

		addDialogComponent(new DialogComponentNumber(max_cycles, "Cycle threshold", 1, 5));
		addDialogComponent(new DialogComponentNumber(gap_open, "Gap open penalty", 0.1, 5));

		addDialogComponent(new DialogComponentNumber(deletion_def_qual,
				"Default quality for deletions (set -1 to disable)", 1, 5));
		addDialogComponent(new DialogComponentNumber(insertion_def_qual,
				"Default quality for insertions (set -1 to disable)", 1, 5));
		addDialogComponent(new DialogComponentNumber(mismatch_def_qual,
				"Default quality for mismatches (set -1 to disable)", 1, 5));
		addDialogComponent(new DialogComponentNumber(indel_context_size, "k-mer context size for indels", 1, 5));
		addDialogComponent(new DialogComponentNumber(mismatch_context_size, "k-mer context size for mismatches", 1, 5));

		addDialogComponent(new DialogComponentNumber(low_qual_tail, "Quality threshold", 1, 5));
		
		addDialogComponent(new DialogComponentOptionalString(m_br_opt_flags, "Optional flags"));

		// disables max value for cycles if cycle covariate is not calculated
//		cycle_cov.addChangeListener(new ChangeListener() {
//
//			public void stateChanged(ChangeEvent e) {
//				max_cycles.setEnabled(cycle_cov.getBooleanValue());
//			}
//		});

		// disables covariate choice
//		create_plots.addChangeListener(new ChangeListener() {
//
//			public void stateChanged(ChangeEvent e) {
//				context_cov.setEnabled(!create_plots.getBooleanValue());
//				cycle_cov.setEnabled(!create_plots.getBooleanValue());
//				rep_len_cov.setEnabled(!create_plots.getBooleanValue());
//				rep_unit_cov.setEnabled(!create_plots.getBooleanValue());
//			}
//		});
	}

	// private void generateProxyOptions(){
	// createNewTab("Proxy options");
	// createNewGroup("General");
	// addDialogComponent(new DialogComponentBoolean(useproxy, "Enable proxy"));
	// addDialogComponent(new DialogComponentString(proxyhost, "Proxy host"));
	// addDialogComponent(new DialogComponentString(proxyport, "Proxy port"));
	// createNewGroup("Authentication");
	// addDialogComponent(new DialogComponentBoolean(useproxyauth, "Enable
	// authentication"));
	// addDialogComponent(new DialogComponentString(proxyuser, "Proxy
	// username"));
	// addDialogComponent(new DialogComponentString(proxypassword, "Proxy
	// password"));
	//
	//
	// useproxy.addChangeListener(new ChangeListener() {
	// public void stateChanged(ChangeEvent e) {
	// if(useproxy.getBooleanValue()){
	// proxyhost.setEnabled(true);
	// proxyport.setEnabled(true);
	// useproxyauth.setEnabled(true);
	// }else{
	// proxyhost.setEnabled(false);
	// proxyport.setEnabled(false);
	// proxyuser.setEnabled(false);
	// proxypassword.setEnabled(false);
	// useproxyauth.setEnabled(false);
	// }
	// }
	// });
	//
	// useproxyauth.addChangeListener(new ChangeListener() {
	// public void stateChanged(ChangeEvent e) {
	// if(useproxy.getBooleanValue() && useproxyauth.getBooleanValue()){
	// proxyuser.setEnabled(true);
	// proxypassword.setEnabled(true);
	// }else{
	// proxypassword.setEnabled(false);
	// proxyuser.setEnabled(false);
	// }
	// }
	// });
}
