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
package de.helmholtz_muenchen.ibis.ngs.samtools;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "SamTools" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more
 * complex dialog please derive directly from
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Sebastian Kopetzky
 * @author Maximilian Hastreiter
 */
public class SamToolsNodeDialog extends HTExecutorNodeDialog {

	protected SamToolsNodeDialog() {
	}

	public void addToolDialogComponents() {
		final SettingsModelString utility 			= new SettingsModelString(SamToolsNodeModel.CFGKEY_UTILITY, "");
		final SettingsModelString samtools 			= new SettingsModelString(SamToolsNodeModel.CFGKEY_SAMTOOLS, "");
		final SettingsModelString refseqfile 		= new SettingsModelString(SamToolsNodeModel.CFGKEY_REFSEQFILE, "");
		final SettingsModelOptionalString Optional 	= new SettingsModelOptionalString(SamToolsNodeModel.CFGKEY_OPTIONAL,"",false);

		
		// calmd
//		final SettingsModelBoolean changeIdentBases = new SettingsModelBoolean(
//				SamToolsNodeModel.CFGKEY_CHANGEIDENTBASES, false);
//		final SettingsModelBoolean useCompression = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_USECOMPRESSION,
//				false);
//		final SettingsModelString compression = new SettingsModelString(SamToolsNodeModel.CFGKEY_COMPRESSION, "");
//		final SettingsModelBoolean inputIsSAM = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_INPUTISSAM, false);
//		final SettingsModelBoolean modifyQual = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MODIFYQUAL, false);
//		final SettingsModelBoolean bqTag = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_BQTAG, false);
//		final SettingsModelBoolean extendedBAQ = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_EXTENDEDBAQ, false);

		// rmdup
		final SettingsModelBoolean removeDup 	= new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_REMOVEDUP, false);
		final SettingsModelBoolean treatPE 		= new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_TREATPE, false);
		// cat
		final SettingsModelBoolean useHeaderSAM = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_USEHEADERSAM,false);
		final SettingsModelString headerSAM 	= new SettingsModelString(SamToolsNodeModel.CFGKEY_HEADERSAM, "");
		final SettingsModelString rehInSAM 		= new SettingsModelString(SamToolsNodeModel.CFGKEY_REHINSAM, "");
		
		// merge
		final SettingsModelBoolean mcompression = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MCOMPRESSION,false);
		final SettingsModelBoolean mforce 		= new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MFORCE, false);
		final SettingsModelBoolean usemhfile 	= new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_USEMHFILE, false);
		final SettingsModelString mhfile 		= new SettingsModelString(SamToolsNodeModel.CFGKEY_MHFILE, "");
		final SettingsModelBoolean msorted 		= new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MSORTED, false);
		final SettingsModelString mregion 		= new SettingsModelString(SamToolsNodeModel.CFGKEY_MREGION, "");
		final SettingsModelBoolean usemregion 	= new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_USEMREGION, false);
		final SettingsModelBoolean mrgtag 		= new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MRGTAG, false);
		final SettingsModelBoolean muncompressed= new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MUNCOMPRESSED,false);
		
		// faidx
		final SettingsModelIntegerBounded blocklength 	= new SettingsModelIntegerBounded(SamToolsNodeModel.CFGKEY_BLOCKLENGTH, 13, 1, Integer.MAX_VALUE);
		final SettingsModelString prefix 				= new SettingsModelString(SamToolsNodeModel.CFGKEY_PREFIX, "phase");
		final SettingsModelIntegerBounded hetphred 		= new SettingsModelIntegerBounded(SamToolsNodeModel.CFGKEY_HETPHRED,37, 1, Integer.MAX_VALUE);
		final SettingsModelIntegerBounded minqual 		= new SettingsModelIntegerBounded(SamToolsNodeModel.CFGKEY_MINQUAL,13, 1, Integer.MAX_VALUE);
		final SettingsModelIntegerBounded maxdepth 		= new SettingsModelIntegerBounded(SamToolsNodeModel.CFGKEY_MAXDEPTH,256, 1, Integer.MAX_VALUE);
		final SettingsModelBoolean fixchimeras 			= new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_FIXCHIMERAS, false);

		addPrefPageSetting(samtools, IBISKNIMENodesPlugin.SAMTOOLS);
		addPrefPageSetting(refseqfile, IBISKNIMENodesPlugin.REF_GENOME);

		createNewGroup("Select utility");
		addDialogComponent(new DialogComponentStringSelection(utility, "Select Utility", "cat", "faidx",
				"fixmate", "flagstat", "idxstats", "merge", "phase", "reheader", "rmdup"));
    	addDialogComponent(new DialogComponentOptionalString(Optional, "Optional Parameters"));
		

		createNewTab("cat");
		setHorizontalPlacement(true);
		createNewGroup("Header SAM file:");
		addDialogComponent(new DialogComponentBoolean(useHeaderSAM, "Use header SAM"));
		addDialogComponent(new DialogComponentFileChooser(headerSAM, "his_id0_samtools", 0, ".sam", ".SAM"));
		setHorizontalPlacement(false);


//		createNewTab("calmd");
//		addDialogComponent(new DialogComponentBoolean(changeIdentBases, "Change ident. bases to '='"));
//		setHorizontalPlacement(true);
//		addDialogComponent(new DialogComponentBoolean(useCompression, "Output as BAM file"));
//		addDialogComponent(new DialogComponentStringSelection(compression, "", "compressed", "uncompressed"));
//		setHorizontalPlacement(false);
//		addDialogComponent(new DialogComponentBoolean(inputIsSAM, "Input is SAM with header"));
//		addDialogComponent(new DialogComponentBoolean(modifyQual, "Modify the quality string"));
//		addDialogComponent(
//				new DialogComponentBoolean(bqTag, "Compute BQ tag if qual. String modif., else cap baseQ by BAQ"));
//		addDialogComponent(
//				new DialogComponentBoolean(extendedBAQ, "Extended BAQ for better sensitivity but lower specificity"));


		createNewTab("merge");
		createNewGroup("Parameters:");
		setHorizontalPlacement(true);
		addDialogComponent(new DialogComponentBoolean(usemhfile, "Use lines of SAM file as header"));
		addDialogComponent(new DialogComponentFileChooser(mhfile, "his_id2_samtools", 0, ".sam", ".SAM"));
		setHorizontalPlacement(false);
		setHorizontalPlacement(true);
		addDialogComponent(new DialogComponentBoolean(usemregion, "Merge files in specified region"));
		addDialogComponent(new DialogComponentString(mregion, ""));
		setHorizontalPlacement(false);
		addDialogComponent(new DialogComponentBoolean(mcompression, "Use zlib lvl 1 for output compression"));
		addDialogComponent(new DialogComponentBoolean(mforce, "Overwrite outputfile if present"));
		addDialogComponent(new DialogComponentBoolean(msorted, "Input alignments are sorted by read names"));
		addDialogComponent(new DialogComponentBoolean(mrgtag, "Attach RG tag to each alignment"));
		addDialogComponent(new DialogComponentBoolean(muncompressed, "Uncompressed BAM output"));

		createNewTab("phase");
		addDialogComponent(new DialogComponentString(prefix, "Prefix of BAM output"));
		addDialogComponent(new DialogComponentNumber(blocklength, "Block length", /* step */ 1));
		addDialogComponent(new DialogComponentNumber(hetphred, "Min. het phred-LOD", /* step */ 1));
		addDialogComponent(new DialogComponentNumber(minqual, "Min. base quality in het calling", /* step */ 1));
		addDialogComponent(new DialogComponentNumber(maxdepth, "Max. read depth", /* step */ 1));
		addDialogComponent(new DialogComponentBoolean(fixchimeras, "Do not attempt to fix chimeric reads"));
		// addDialogComponent(new DialogComponentBoolean(dropambig, "Drop reads
		// with ambiguous phase"));

		createNewTab("reheader");
		createNewGroup("Header SAM file:");
		addDialogComponent(new DialogComponentFileChooser(rehInSAM, "Header SAM file", 0, ".sam", ".SAM"));
		// createNewGroup("BAM file:");
		// addDialogComponent(new DialogComponentFileChooser(rehInBAM, "BAM
		// file", 0,".bam", ".BAM"));

		createNewTab("rmdup");
		addDialogComponent(new DialogComponentBoolean(removeDup, "Remove duplicate for single-end reads"));
		addDialogComponent(new DialogComponentBoolean(treatPE, "Treat paired-end reads and single-end reads."));

		/** Checkboxes for the optional arguments **/

		// checkboxes for calmd tab
		utility.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				setEnabled(true, "cat");
				
//				// calmd
//				if (utility.getStringValue().equals("calmd")) {
//					setEnabled(true, "calmd");
//				} else {
//					setEnabled(false, "calmd");
//
//				}
				// rmdup
				if (utility.getStringValue().equals("rmdup")) {
					setEnabled(true, "rmdup");
				} else {
					setEnabled(false, "rmdup");
				}
				// cat
				if (utility.getStringValue().equals("cat")) {
					setEnabled(true, "cat");
				} else {
					setEnabled(false, "cat");
				}
				// reheader
				if (utility.getStringValue().equals("reheader")) {
					setEnabled(true, "reheader");
				} else {
					setEnabled(false, "reheader");
				}
				// merge
				if (utility.getStringValue().equals("merge")) {
					setEnabled(true, "merge");
				} else {
					setEnabled(false, "merge");
				}
				
				// phase
				if (utility.getStringValue().equals("phase")) {
					setEnabled(true, "phase");

				} else {
					setEnabled(false, "phase");
				}

			}
		});

//		useCompression.addChangeListener(new ChangeListener() {
//			public void stateChanged(ChangeEvent e) {
//				if (useCompression.getBooleanValue()) {
//					compression.setEnabled(true);
//				} else {
//					compression.setEnabled(false);
//				}
//			}
//		});
		useHeaderSAM.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if (useHeaderSAM.getBooleanValue()) {
					headerSAM.setEnabled(true);
				} else {
					headerSAM.setEnabled(false);
				}
			}
		});
		usemhfile.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if (usemhfile.getBooleanValue()) {
					mhfile.setEnabled(true);
				} else {
					mhfile.setEnabled(false);
				}
			}
		});
		usemregion.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if (usemregion.getBooleanValue()) {
					mregion.setEnabled(true);
				} else {
					mregion.setEnabled(false);
				}
			}
		});

	}
}
