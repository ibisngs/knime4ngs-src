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

package de.helmholtz_muenchen.ibis.ngs.bowtie2;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentNumberEdit;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FastQCell;

/**
 * <code>NodeDialog</code> for the "Bowtie2" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Jan Quell 
 */
public class Bowtie2NodeDialog extends HTExecutorNodeDialog {

    /**
     * New pane for configuring the Bowtie2 node.
     */
    protected Bowtie2NodeDialog() {
    	super(FastQCell.TYPE.getPreferredValueClass(),FastQCell.TYPE.getPreferredValueClass());
    }
    	
    public void addToolDialogComponents() {
        	
    	final SettingsModelString installpath = new SettingsModelString(Bowtie2NodeModel.CFGKEY_INSTALLPATH,"");
    	final SettingsModelString refseqfile = new SettingsModelString(Bowtie2NodeModel.CFGKEY_REFSEQFILE,"");
    	final SettingsModelBoolean noauto = new SettingsModelBoolean(Bowtie2NodeModel.CFGKEY_NOAUTO, true);
    	final SettingsModelBoolean packed = new SettingsModelBoolean(Bowtie2NodeModel.CFGKEY_PACKED, false);
    	final SettingsModelIntegerBounded bmax = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_BMAX,Bowtie2NodeModel.DEFAULT_BMAX,0,Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded dcv = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_DCV,Bowtie2NodeModel.DEFAULT_DCV,1,4096);
    	final SettingsModelBoolean nodc = new SettingsModelBoolean(Bowtie2NodeModel.CFGKEY_NODC, false);
    	final SettingsModelIntegerBounded offrate = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_OFFRATE,Bowtie2NodeModel.DEFAULT_OFFRATE,0,Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded ftabchars = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_FTABCHARS,Bowtie2NodeModel.DEFAULT_FTABCHARS,1,Integer.MAX_VALUE);
    	final SettingsModelBoolean useskip = new SettingsModelBoolean(Bowtie2NodeModel.CFGKEY_USESKIP, false);
    	final SettingsModelIntegerBounded skip = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_SKIP,Bowtie2NodeModel.DEFAULT_SKIP,0,Integer.MAX_VALUE);
    	final SettingsModelBoolean useupto = new SettingsModelBoolean(Bowtie2NodeModel.CFGKEY_USEUPTO, false);
    	final SettingsModelIntegerBounded upto = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_UPTO,Bowtie2NodeModel.DEFAULT_UPTO,0,Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded trim5 = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_TRIM5,Bowtie2NodeModel.DEFAULT_TRIM5,0,Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded trim3 = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_TRIM3,Bowtie2NodeModel.DEFAULT_TRIM3,0,Integer.MAX_VALUE);
    	final SettingsModelString quals = new SettingsModelString(Bowtie2NodeModel.CFGKEY_QUALS,"");
    	final SettingsModelBoolean usepreset = new SettingsModelBoolean(Bowtie2NodeModel.CFGKEY_USEPRESET, true);
    	final SettingsModelString preset = new SettingsModelString(Bowtie2NodeModel.CFGKEY_PRESET,"");
    	final SettingsModelIntegerBounded n = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_N,Bowtie2NodeModel.DEFAULT_N,0,1);
    	final SettingsModelIntegerBounded l = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_L,Bowtie2NodeModel.DEFAULT_L,4,31);
    	final SettingsModelDoubleBounded i1 = new SettingsModelDoubleBounded(Bowtie2NodeModel.CFGKEY_I1,Bowtie2NodeModel.DEFAULT_I1,0,Double.MAX_VALUE);
    	final SettingsModelDoubleBounded i2 = new SettingsModelDoubleBounded(Bowtie2NodeModel.CFGKEY_I2,Bowtie2NodeModel.DEFAULT_I2,0,Double.MAX_VALUE);
    	final SettingsModelDoubleBounded nceil1 = new SettingsModelDoubleBounded(Bowtie2NodeModel.CFGKEY_NCEIL1,Bowtie2NodeModel.DEFAULT_NCEIL1,0,Double.MAX_VALUE);
    	final SettingsModelDoubleBounded nceil2 = new SettingsModelDoubleBounded(Bowtie2NodeModel.CFGKEY_NCEIL2,Bowtie2NodeModel.DEFAULT_NCEIL2,0,Double.MAX_VALUE);
    	final SettingsModelIntegerBounded dpad = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_DPAD,Bowtie2NodeModel.DEFAULT_DPAD,0,Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded gbar = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_GBAR,Bowtie2NodeModel.DEFAULT_GBAR,0,Integer.MAX_VALUE);
    	final SettingsModelBoolean ignorequals = new SettingsModelBoolean(Bowtie2NodeModel.CFGKEY_IGNOREQUALS, false);
    	final SettingsModelBoolean nofw = new SettingsModelBoolean(Bowtie2NodeModel.CFGKEY_NOFW, false);
    	final SettingsModelBoolean norc = new SettingsModelBoolean(Bowtie2NodeModel.CFGKEY_NORC, false);
    	final SettingsModelString alignmenttype = new SettingsModelString(Bowtie2NodeModel.CFGKEY_ALIGNMENTTYPE,"");
    	final SettingsModelIntegerBounded ma = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_MA,Bowtie2NodeModel.DEFAULT_MA,0,Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded mp = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_MP,Bowtie2NodeModel.DEFAULT_MP,0,Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded np = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_NP,Bowtie2NodeModel.DEFAULT_NP,0,Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded rdg1 = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_RDG1,Bowtie2NodeModel.DEFAULT_RDG1,0,Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded rdg2 = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_RDG2,Bowtie2NodeModel.DEFAULT_RDG2,0,Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded rfg1 = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_RFG1,Bowtie2NodeModel.DEFAULT_RFG1,0,Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded rfg2 = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_RFG2,Bowtie2NodeModel.DEFAULT_RFG2,0,Integer.MAX_VALUE);
    	final SettingsModelDoubleBounded scoremin1 = new SettingsModelDoubleBounded(Bowtie2NodeModel.CFGKEY_SCOREMIN1,Bowtie2NodeModel.DEFAULT_SCOREMIN1,-Double.MAX_VALUE,Double.MAX_VALUE);
    	final SettingsModelDoubleBounded scoremin2 = new SettingsModelDoubleBounded(Bowtie2NodeModel.CFGKEY_SCOREMIN2,Bowtie2NodeModel.DEFAULT_SCOREMIN2,-Double.MAX_VALUE,Double.MAX_VALUE);
    	final SettingsModelString reporting1 = new SettingsModelString(Bowtie2NodeModel.CFGKEY_REPORTING1,"");
    	final SettingsModelIntegerBounded reporting2 = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_REPORTING2,Bowtie2NodeModel.DEFAULT_REPORTING2,0,Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded d = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_D,Bowtie2NodeModel.DEFAULT_D,0,Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded r = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_R,Bowtie2NodeModel.DEFAULT_R,0,Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded minins = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_MININS,Bowtie2NodeModel.DEFAULT_MININS,0,Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded maxins = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_MAXINS,Bowtie2NodeModel.DEFAULT_MAXINS,0,Integer.MAX_VALUE);
    	final SettingsModelString ff = new SettingsModelString(Bowtie2NodeModel.CFGKEY_FF,"");
    	final SettingsModelBoolean nomixed = new SettingsModelBoolean(Bowtie2NodeModel.CFGKEY_NOMIXED, false);
    	final SettingsModelBoolean nodiscordant = new SettingsModelBoolean(Bowtie2NodeModel.CFGKEY_NODISCORDANT, false);
    	final SettingsModelBoolean nodovetail = new SettingsModelBoolean(Bowtie2NodeModel.CFGKEY_NODOVETAIL, false);
    	final SettingsModelBoolean nocontain = new SettingsModelBoolean(Bowtie2NodeModel.CFGKEY_NOCONTAIN, false);
    	final SettingsModelBoolean nooverlap = new SettingsModelBoolean(Bowtie2NodeModel.CFGKEY_NOOVERLAP, false);
    	final SettingsModelIntegerBounded threads = new SettingsModelIntegerBounded(Bowtie2NodeModel.CFGKEY_THREADS,Bowtie2NodeModel.DEFAULT_THREADS,1,Integer.MAX_VALUE);
    	final SettingsModelBoolean recorder = new SettingsModelBoolean(Bowtie2NodeModel.CFGKEY_RECORDER, false);
    	final SettingsModelBoolean mm = new SettingsModelBoolean(Bowtie2NodeModel.CFGKEY_MM, true);
    	final SettingsModelBoolean qcfilter = new SettingsModelBoolean(Bowtie2NodeModel.CFGKEY_QCFILTER, false);
    	
    	
    	addPrefPageSetting(installpath, IBISKNIMENodesPlugin.BOWTIE2);
    	addPrefPageSetting(refseqfile, IBISKNIMENodesPlugin.REF_GENOME);
    	
    	
    	packed.setEnabled(false);
    	bmax.setEnabled(false);
    	dcv.setEnabled(false);
    	skip.setEnabled(false);
    	upto.setEnabled(false);
    	n.setEnabled(false);
    	l.setEnabled(false);
    	i1.setEnabled(false);
    	i2.setEnabled(false);
    	preset.setStringValue("sensitive");
    	alignmenttype.setStringValue("entire read must align (no clipping)");
    	ma.setEnabled(false);
    	reporting2.setEnabled(false);
    	d.setEnabled(false);
    	r.setEnabled(false);
    	minins.setEnabled(false);
    	maxins.setEnabled(false);
    	ff.setEnabled(false);
    	nomixed.setEnabled(false);
    	nodiscordant.setEnabled(false);
    	nodovetail.setEnabled(false);
    	nocontain.setEnabled(false);
    	nooverlap.setEnabled(false);
    	
//    	createNewGroup("Bowtie2 binary");
//    	addDialogComponent(new DialogComponentFileChooser(installpath, "his_bow_id", 0, ""));
//    	createNewGroup("Reference sequence: FastA file (e.g. genome)");
//    	addDialogComponent(new DialogComponentFileChooser(refseqfile, "his1_bow_id", 0, ""));
    	createNewGroup("Indexing parameters");
    	addDialogComponent(new DialogComponentBoolean(noauto, "Automatically select value for parameters according to available memory."));
    	addDialogComponent(new DialogComponentBoolean(packed, "Use a packed (2-bits-per-nucleotide) representation for DNA strings."));
    	addDialogComponent(new DialogComponentNumber(bmax,"The maximum number of suffixes allowed in a block:", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumberEdit(dcv, "Use <int> as the period for the difference-cover sample (must be a power of 2 no greater than 4096):"));
    	addDialogComponent(new DialogComponentBoolean(nodc, "Disable use of the difference-cover sample (algorithm becomes quadratic)"));
    	addDialogComponent(new DialogComponentNumber(offrate,"Mark every 2^<int> rows. Marking more rows makes lookups faster, but requires more memory:", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(ftabchars,"Use the first <int> characters of the query to calculate an initial Burrows-Wheeler:", /*step*/ 1));
    	createNewTab("Alignment Parameters");
    	createNewGroup("Input");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(useskip, "Skip the first <int> reads/pairs in the input:"));
    	addDialogComponent(new DialogComponentNumber(skip,"", /*step*/ 1));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(useupto, "Stop after first <int> reads/pairs:"));
    	addDialogComponent(new DialogComponentNumber(upto,"", /*step*/ 1));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentNumber(trim5,
    			"Trim <int> bases from 5'/left end of reads:", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(trim3,
    			"Trim <int> bases from 3'/right end of reads:", /*step*/ 1));
    	addDialogComponent(new DialogComponentStringSelection(quals,"Select quality score type:", 
    			"Phred+33", "Phred+64", "Convert from Solexa to Phred", "Encoded as space-delimited integers"));
    	createNewGroup("Presets");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usepreset, "Use a set of preselected parameters:"));
    	addDialogComponent(new DialogComponentStringSelection(preset, "", "very-fast", "fast", "sensitive", "very-sensitive"));
    	setHorizontalPlacement(false);
    	createNewGroup("Alignment");
    	addDialogComponent(new DialogComponentNumber(n,"Max <int> mismatches in seed alignment (can be 0 or 1):", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(l,"Length of seed substrings (must be >3, <32):", /*step*/ 1));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentNumber(i1,"Interval between seed substrings w/r/t read len:", /*step*/ 0.01));
    	addDialogComponent(new DialogComponentNumber(i2," - ", /*step*/ 0.01));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentNumber(nceil1,"Func for max # non-A/C/G/Ts permitted in aln:", /*step*/ 0.01));
    	addDialogComponent(new DialogComponentNumber(nceil2," - ", /*step*/ 0.01));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentNumber(dpad,
    			"Include <int> extra ref chars on sides of DP table:", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(gbar,
    			"Disallow gaps within <int> nucs of read extremes:", /*step*/ 1));
    	addDialogComponent(new DialogComponentBoolean(ignorequals, "Treat all quality values as 30 on Phred scale."));
    	addDialogComponent(new DialogComponentBoolean(nofw, "Do not align forward (original) version of read."));
    	addDialogComponent(new DialogComponentBoolean(norc, "Do not align reverse-complement version of read."));
    	addDialogComponent(new DialogComponentStringSelection(alignmenttype,"Select an alignment type:",
    			"entire read must align (no clipping)", "local alignment (ends might be soft clipped)"));
    	createNewTab("Further parameters");
    	createNewGroup("Scoring");
    	addDialogComponent(new DialogComponentNumber(ma, "Set the match bonus (<int> is added to the alignment score for each position):", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(mp,
    			"Max penalty for mismatch (lower qual = lower penalty):", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(np,
    			"Penalty for non-A/C/G/Ts in read/ref:", /*step*/ 1));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentNumber(rdg1,
    			"Read gap open penalty:", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(rdg2,
    			"Read gap extend penalty:", /*step*/ 1));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentNumber(rfg1,
    			"Reference gap open penalty:", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(rfg2,
    			"Reference gap extend penalty:", /*step*/ 1));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentNumber(scoremin1,"Min acceptable alignment score w/r/t read length:", /*step*/ 0.1));
    	addDialogComponent(new DialogComponentNumber(scoremin2," - ", /*step*/ 0.1));
    	setHorizontalPlacement(false);
    	createNewGroup("Reporting");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentStringSelection(reporting1,"Choose an alignment type:","Look for multiple alignments, report best, with MAPQ", 
    			"Report up to <int> alns per read; MAPQ not meaningful:", "Report all alignments; very slow, MAPQ not meaningful"));
    	addDialogComponent(new DialogComponentNumber(reporting2, "", /*step*/ 1));
    	setHorizontalPlacement(false);
    	createNewGroup("Effort");
    	addDialogComponent(new DialogComponentNumber(d, "Give up extending after <int> failed extends in a row:", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(r, "For reads with repetitive seeds, try <int> sets of seeds:", /*step*/ 1));
    	createNewGroup("Paired-end");
    	addDialogComponent(new DialogComponentNumber(minins, "Minimum fragment length:", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(maxins, "Maximum fragment length:", /*step*/ 1));
    	addDialogComponent(new DialogComponentStringSelection(ff,"Select upstream/downstream mate orientations in the alignment:","forward/reverse", 
    			"reverse/forward", "forward/forward"));
    	addDialogComponent(new DialogComponentBoolean(nomixed, "Suppress unpaired alignments for paired reads."));
    	addDialogComponent(new DialogComponentBoolean(nodiscordant, "Suppress discordant alignments for paired reads."));
    	addDialogComponent(new DialogComponentBoolean(nodovetail, "Not concordant when mates extend past each other."));
    	addDialogComponent(new DialogComponentBoolean(nocontain, "Not concordant when one mate alignment contains other."));
    	addDialogComponent(new DialogComponentBoolean(nooverlap, "Not concordant when mates overlap at all."));
    	createNewGroup("Performance");
    	addDialogComponent(new DialogComponentNumber(threads, "Number of alignment threads to launch (cores to use):", /*step*/ 1));
    	addDialogComponent(new DialogComponentBoolean(recorder, "Force SAM output order to match order of input reads."));
    	addDialogComponent(new DialogComponentBoolean(mm, "Use memory-mapped I/O for index (many threads of bowtie can share)."));
    	createNewGroup("Other");
    	addDialogComponent(new DialogComponentBoolean(qcfilter, "Filter out reads that are bad according to QSEQ filter."));
    	
    	
    	noauto.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				packed.setEnabled(!noauto.getBooleanValue());
				bmax.setEnabled(!noauto.getBooleanValue());
				dcv.setEnabled(!noauto.getBooleanValue());
			}
		});
    	

    	
    	useskip.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				skip.setEnabled(useskip.getBooleanValue());
			}
		});
    	
    	useupto.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				upto.setEnabled(useupto.getBooleanValue());
			}
		});
    	
    	usepreset.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				preset.setEnabled(usepreset.getBooleanValue());
				d.setEnabled(!usepreset.getBooleanValue());
				r.setEnabled(!usepreset.getBooleanValue());
				n.setEnabled(!usepreset.getBooleanValue());
				l.setEnabled(!usepreset.getBooleanValue());
				i1.setEnabled(!usepreset.getBooleanValue());
				i2.setEnabled(!usepreset.getBooleanValue());
			}
		});
    	
    	alignmenttype.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				ma.setEnabled(false);
				scoremin1.setDoubleValue(-0.6);
				scoremin2.setDoubleValue(-0.6);
				if(alignmenttype.getStringValue().equals("local alignment (ends might be soft clipped)")) {
					ma.setEnabled(true);
					scoremin1.setDoubleValue(20);
					scoremin2.setDoubleValue(8);
				}
			}
		});
    	
    	reporting1.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				reporting2.setEnabled(reporting1.getStringValue().equals("Report up to <int> alns per read; MAPQ not meaningful:"));
			}
		});
    	
//    	readType.addChangeListener(new ChangeListener() {
//			public void stateChanged(ChangeEvent e) {
//				minins.setEnabled(readType.getStringValue().equals("paired-end"));
//		    	maxins.setEnabled(readType.getStringValue().equals("paired-end"));
//		    	ff.setEnabled(readType.getStringValue().equals("paired-end"));
//		    	nomixed.setEnabled(readType.getStringValue().equals("paired-end"));
//		    	nodiscordant.setEnabled(readType.getStringValue().equals("paired-end"));
//		    	nodovetail.setEnabled(readType.getStringValue().equals("paired-end"));
//		    	nocontain.setEnabled(readType.getStringValue().equals("paired-end"));
//		    	nooverlap.setEnabled(readType.getStringValue().equals("paired-end"));
//			}
//		});
    	
    }
}

