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
package de.helmholtz_muenchen.ibis.ngs.gatkrealignment;

import javax.swing.JFileChooser;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "GATKRealignment" Node.
 * 
 * 
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more
 * complex dialog please derive directly from
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Marie-Sophie Friedl
 * @author Maximilian Hastreiter
 * @author Tim Jeske
 */
public class GATKRealignmentNodeDialog extends HTExecutorNodeDialog {

	protected GATKRealignmentNodeDialog() {
	};

	public void addToolDialogComponents() {

		//pref page
		final SettingsModelString gatk = new SettingsModelString(GATKRealignmentNodeModel.CFGKEY_GATK, "");
		final SettingsModelString reffile = new SettingsModelString(GATKRealignmentNodeModel.CFGKEY_REF_GENOME, "");
		final SettingsModelString phase1_1000G_file = new SettingsModelString(
				GATKRealignmentNodeModel.CFGKEY_PHASE1_1000G_FILE, GATKRealignmentNodeModel.DEF_PHASE1_1000G_FILE);
		final SettingsModelString mills_1000G_file = new SettingsModelString(
				GATKRealignmentNodeModel.CFGKEY_MILLS_1000G_FILE, GATKRealignmentNodeModel.DEF_MILLS_1000G_FILE);
		
		final SettingsModelBoolean use_phase1_1000G = new SettingsModelBoolean(
				GATKRealignmentNodeModel.CFGKEY_USE_PHASE1_1000G, GATKRealignmentNodeModel.DEF_USE_PHASE1_1000G);
		final SettingsModelBoolean use_mills_1000G = new SettingsModelBoolean(
				GATKRealignmentNodeModel.CFGKEY_USE_MILLS_1000G, GATKRealignmentNodeModel.DEF_USE_MILLS_1000G);
		final SettingsModelBoolean use_interval = new SettingsModelBoolean(GATKRealignmentNodeModel.CFGKEY_USE_INTERVAL,
				GATKRealignmentNodeModel.DEF_USE_INTERVAL);
		final SettingsModelString interval_file = new SettingsModelString(GATKRealignmentNodeModel.CFGKEY_INTERVAL_FILE,
				GATKRealignmentNodeModel.DEF_INTERVAL_FILE);
		final SettingsModelIntegerBounded num_threads = new SettingsModelIntegerBounded(
				GATKRealignmentNodeModel.CFGKEY_NUM_THREADS, GATKRealignmentNodeModel.DEF_NUM_THREADS,
				GATKRealignmentNodeModel.MIN_NUM_THREADS, GATKRealignmentNodeModel.MAX_NUM_THREADS);
		final SettingsModelIntegerBounded memory_usage = new SettingsModelIntegerBounded(
				GATKRealignmentNodeModel.CFGKEY_JAVAMEMORY, GATKRealignmentNodeModel.DEF_NUM_JAVAMEMORY,
				GATKRealignmentNodeModel.MIN_NUM_JAVAMEMORY, GATKRealignmentNodeModel.MAX_NUM_JAVAMEMORY);

		final SettingsModelIntegerBounded max_interval = new SettingsModelIntegerBounded(
				GATKRealignmentNodeModel.CFGKEY_MAX_INTERVAL, GATKRealignmentNodeModel.DEF_MAX_INTERVAL,
				GATKRealignmentNodeModel.MIN_MAX_INTERVAL, GATKRealignmentNodeModel.MAX_MAX_INTERVAL);
		final SettingsModelIntegerBounded min_reads = new SettingsModelIntegerBounded(
				GATKRealignmentNodeModel.CFGKEY_MIN_READS, GATKRealignmentNodeModel.DEF_MIN_READS,
				GATKRealignmentNodeModel.MIN_MIN_READS, GATKRealignmentNodeModel.MAX_MIN_READS);
		final SettingsModelIntegerBounded window = new SettingsModelIntegerBounded(
				GATKRealignmentNodeModel.CFGKEY_WINDOW, GATKRealignmentNodeModel.DEF_WINDOW,
				GATKRealignmentNodeModel.MIN_WINDOW, GATKRealignmentNodeModel.MAX_WINDOW);
		final SettingsModelDoubleBounded mismatch = new SettingsModelDoubleBounded(
				GATKRealignmentNodeModel.CFGKEY_MISMATCH, GATKRealignmentNodeModel.DEF_MISMATCH,
				GATKRealignmentNodeModel.MIN_MISMATCH, GATKRealignmentNodeModel.MAX_MISMATCH);
		final SettingsModelOptionalString m_tc_opt_flags = new SettingsModelOptionalString(
				GATKRealignmentNodeModel.CFGKEY_TC_OPT_FLAGS, "", false);

		final SettingsModelString consensus_model = new SettingsModelString(
				GATKRealignmentNodeModel.CFGKEY_CONSENSUS_MODEL, GATKRealignmentNodeModel.DEF_CONSENSUS_MODEL);
		final SettingsModelDoubleBounded lod_threshold = new SettingsModelDoubleBounded(
				GATKRealignmentNodeModel.CFGKEY_LOD_THRESHOLD, GATKRealignmentNodeModel.DEF_LOD_THRESHOLD,
				GATKRealignmentNodeModel.MIN_LOD_THRESHOLD, GATKRealignmentNodeModel.MAX_LOD_THRESHOLD);
		final SettingsModelDoubleBounded entropy = new SettingsModelDoubleBounded(
				GATKRealignmentNodeModel.CFGKEY_ENTROPY, GATKRealignmentNodeModel.DEF_ENTROPY,
				GATKRealignmentNodeModel.MIN_ENTROPY, GATKRealignmentNodeModel.MAX_ENTROPY);
		final SettingsModelIntegerBounded max_consensuses = new SettingsModelIntegerBounded(
				GATKRealignmentNodeModel.CFGKEY_MAX_CONSENSUSES, GATKRealignmentNodeModel.DEF_MAX_CONSENSUSES,
				GATKRealignmentNodeModel.MIN_MAX_CONSENSUSES, GATKRealignmentNodeModel.MAX_MAX_CONSENSUSES);
		final SettingsModelIntegerBounded max_isize = new SettingsModelIntegerBounded(
				GATKRealignmentNodeModel.CFGKEY_MAX_ISIZE, GATKRealignmentNodeModel.DEF_MAX_ISIZE,
				GATKRealignmentNodeModel.MIN_MAX_ISIZE, GATKRealignmentNodeModel.MAX_MAX_ISIZE);
		final SettingsModelIntegerBounded max_pos_move = new SettingsModelIntegerBounded(
				GATKRealignmentNodeModel.CFGKEY_MAX_POS_MOVE, GATKRealignmentNodeModel.DEF_MAX_POS_MOVE,
				GATKRealignmentNodeModel.MIN_MAX_POS_MOVE, GATKRealignmentNodeModel.MAX_MAX_POS_MOVE);
		final SettingsModelIntegerBounded max_reads_cons = new SettingsModelIntegerBounded(
				GATKRealignmentNodeModel.CFGKEY_MAX_READS_CONS, GATKRealignmentNodeModel.DEF_MAX_READS_CONS,
				GATKRealignmentNodeModel.MIN_MAX_READS_CONS, GATKRealignmentNodeModel.MAX_MAX_READS_CONS);
		final SettingsModelIntegerBounded max_reads_realign = new SettingsModelIntegerBounded(
				GATKRealignmentNodeModel.CFGKEY_MAX_READS_REALIGN, GATKRealignmentNodeModel.DEF_MAX_READS_REALIGN,
				GATKRealignmentNodeModel.MIN_MAX_READS_REALIGN, GATKRealignmentNodeModel.MAX_MAX_READS_REALIGN);
		final SettingsModelBoolean alignment_tag = new SettingsModelBoolean(
				GATKRealignmentNodeModel.CFGKEY_ALIGNMENT_TAG, GATKRealignmentNodeModel.DEF_ALIGNMENT_TAG);
		final SettingsModelOptionalString m_ir_opt_flags = new SettingsModelOptionalString(
				GATKRealignmentNodeModel.CFGKEY_IR_OPT_FLAGS, "", false);
		
		// Proxy options
		// private final SettingsModelBoolean useproxy = new
		// SettingsModelBoolean(
		// GATKRealignmentNodeModel.CFGKEY_USEPROXY, false);
		// final SettingsModelString proxyhost = new SettingsModelString(
		// GATKRealignmentNodeModel.CFGKEY_PROXYHOST, null);
		// final SettingsModelString proxyport = new SettingsModelString(
		// GATKRealignmentNodeModel.CFGKEY_PROXYPORT, null);
		// private final SettingsModelBoolean useproxyauth = new
		// SettingsModelBoolean(
		// GATKRealignmentNodeModel.CFGKEY_USEPROXYAUTH, false);
		// final SettingsModelString proxyuser = new SettingsModelString(
		// GATKRealignmentNodeModel.CFGKEY_PROXYUSER, null);
		// final SettingsModelString proxypassword = new SettingsModelString(
		// GATKRealignmentNodeModel.CFGKEY_PROXYPASSWORD, null);

		addPrefPageSetting(gatk, IBISKNIMENodesPlugin.GATK);
		addPrefPageSetting(reffile, IBISKNIMENodesPlugin.REF_GENOME);
		addPrefPageSetting(phase1_1000G_file, IBISKNIMENodesPlugin.RES_1000G_INDELS);
		addPrefPageSetting(mills_1000G_file, IBISKNIMENodesPlugin.RES_MILLS);

		// sets of known indels from database for realignment
		createNewGroup("Sets of known indels");
		addDialogComponent(new DialogComponentBoolean(use_phase1_1000G,
				"Use 1000 genomes phase1 indel set?"));

		addDialogComponent(new DialogComponentBoolean(use_mills_1000G,
				"Use Mills and 1000 genomes gold standard indel set?"));

		// interval file for realignment
		createNewGroup("Interval for realignment");
		addDialogComponent(new DialogComponentBoolean(use_interval, "Restrict realignment to certain genomic regions?"));
		interval_file.setEnabled(false);
		addDialogComponent(new DialogComponentFileChooser(interval_file, "ifile", JFileChooser.OPEN_DIALOG, false,
				".bed", ".intervals"));

		// add change listener that enables file chooser for interval file
		use_interval.addChangeListener(new ChangeListener() {

			public void stateChanged(ChangeEvent e) {
				interval_file.setEnabled(use_interval.getBooleanValue());
			}
		});

		createNewGroup("General Options");
		addDialogComponent(new DialogComponentNumber(num_threads, "Threads", 1, 5));
		addDialogComponent(new DialogComponentNumber(memory_usage, "Java Memory (GB) per thread", 1, 5));

		//TargetCreator tab
		createNewTab("RealignerTargetCreator");
		addDialogComponent(new DialogComponentNumber(max_interval, "Maximum length of realignment interval", 1, 5));
		addDialogComponent(new DialogComponentNumber(min_reads, "Minimum number of reads for entropy calculation", 1, 5));
		addDialogComponent(new DialogComponentNumber(mismatch,
				"Fraction of mismatching base qualities", 0.01, 5));
		addDialogComponent(new DialogComponentNumber(window, "Window size for clustering SNPs", 1, 5));
		addDialogComponent(new DialogComponentOptionalString(m_tc_opt_flags, "Optional flags"));
		
		//IndelRealigner tab
		createNewTab("IndelRealigner");
		addDialogComponent(new DialogComponentStringSelection(consensus_model, "Consensus determination model",
				GATKRealignmentNodeModel.VALUES_CONSENSUS_MODEL));
		addDialogComponent(new DialogComponentNumber(lod_threshold, "LOD Threshold", 0.1, 5));
		addDialogComponent(new DialogComponentNumber(entropy, "Entropy threshold", 0.01, 5));
		addDialogComponent(new DialogComponentNumber(max_consensuses, "Consensus threshold", 1, 5));
		addDialogComponent(new DialogComponentNumber(max_isize, "Insert size threshold", 10, 5));
		addDialogComponent(new DialogComponentNumber(max_pos_move, "Read shift threshold", 1, 5));
		addDialogComponent(new DialogComponentNumber(max_reads_cons, "Reads for consensus threshold", 1, 5));
		addDialogComponent(new DialogComponentNumber(max_reads_realign, "Reads for realignemnt threshold", 100, 5));
		addDialogComponent(new DialogComponentBoolean(alignment_tag, "Do not output original cigar string?"));
		addDialogComponent(new DialogComponentOptionalString(m_ir_opt_flags, "Optional flags"));
	}

	// private void createProxyOptions() {
	// createNewTab("Proxy options");
	// createNewGroup("General");
	// addDialogComponent(new DialogComponentBoolean(useproxy, "Enable proxy"));
	// addDialogComponent(new DialogComponentString(proxyhost, "Proxy host"));
	// addDialogComponent(new DialogComponentString(proxyport, "Proxy port"));
	// createNewGroup("Authentication");
	// addDialogComponent(new DialogComponentBoolean(useproxyauth,
	// "Enable authentication"));
	// addDialogComponent(new DialogComponentString(proxyuser,
	// "Proxy username"));
	// addDialogComponent(new DialogComponentString(proxypassword,
	// "Proxy password"));
	//
	// useproxy.addChangeListener(new ChangeListener() {
	// public void stateChanged(ChangeEvent e) {
	// if (useproxy.getBooleanValue()) {
	// proxyhost.setEnabled(true);
	// proxyport.setEnabled(true);
	// useproxyauth.setEnabled(true);
	// } else {
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
	// if (useproxy.getBooleanValue()
	// && useproxyauth.getBooleanValue()) {
	// proxyuser.setEnabled(true);
	// proxypassword.setEnabled(true);
	// } else {
	// proxypassword.setEnabled(false);
	// proxyuser.setEnabled(false);
	// }
	// }
	// });
}