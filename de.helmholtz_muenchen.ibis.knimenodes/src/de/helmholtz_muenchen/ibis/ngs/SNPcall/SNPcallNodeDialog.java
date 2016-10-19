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

package de.helmholtz_muenchen.ibis.ngs.SNPcall;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "SNPcall" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Maximilian Hastreiter
 * @author Sebastian Kopetzky
 * @author Jan Quell
 * 
 */
public class SNPcallNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring SNPcall node dialog.
     * This is just a suggestion to demonstrate possible default dialog
     * components.
     */
    protected SNPcallNodeDialog() {
        
    	/* 
    	 * Input Options:
       -6           assume the quality is in the Illumina-1.3+ encoding
       -A           count anomalous read pairs
       -B           disable BAQ computation
       -b FILE      list of input BAM files [null]
       -C INT       parameter for adjusting mapQ; 0 to disable [0]
       -d INT       max per-BAM depth to avoid excessive memory usage [250]
       -E           extended BAQ for higher sensitivity but lower specificity
       -f FILE      faidx indexed reference sequence file [null]
       -G FILE      exclude read groups listed in FILE [null]
       -l FILE      list of positions (chr pos) or regions (BED) [null]
       -M INT       cap mapping quality at INT [60]
       -r STR       region in which pileup is generated [null]
       -R           ignore RG tags
       -q INT       skip alignments with mapQ smaller than INT [0]
       -Q INT       skip bases with baseQ/BAQ smaller than INT [13]

##SNP Call specific parameters
       -e INT       Phred-scaled gap extension seq error probability [20]
       -F FLOAT     minimum fraction of gapped reads for candidates [0.002]
       -h INT       coefficient for homopolymer errors [100]
       -I           do not perform indel calling
       -L INT       max per-sample depth for INDEL calling [250]
       -m INT       minimum gapped reads for indel candidates [1]
       -o INT       Phred-scaled gap open sequencing error probability [40]
       -P STR       comma separated list of platforms for indels [all]
              
 ##varFilter
Options: -Q INT    minimum RMS mapping quality for SNPs [10]
         -d INT    minimum read depth [2]
         -D INT    maximum read depth [10000000]
         -a INT    minimum number of alternate bases [2]
         -w INT    SNP within INT bp around a gap to be filtered [3]
         -W INT    window size for filtering adjacent gaps [10]
         -1 FLOAT  min P-value for strand bias (given PV4) [0.0001]
         -2 FLOAT  min P-value for baseQ bias [1e-100]
         -3 FLOAT  min P-value for mapQ bias [0]
         -4 FLOAT  min P-value for end distance bias [0.0001]
         -e FLOAT  min P-value for HWE (plus F<0) [0.0001]
         -p        print filtered variants
    	 */
    	
    	final SettingsModelString bedfile=new SettingsModelString(SNPcallNodeModel.CFGKEY_BEDFILE,null);
		final SettingsModelBoolean ifbedfile = new SettingsModelBoolean(SNPcallNodeModel.CFGKEY_IFBEDFILE, false);
    	bedfile.setEnabled(false);
    	final SettingsModelString excludereadfile=new SettingsModelString(SNPcallNodeModel.CFGKEY_EXCLUDEREADS,null);
		final SettingsModelBoolean ifreadfile = new SettingsModelBoolean(SNPcallNodeModel.CFGKEY_IFREADFILE, false);
    	excludereadfile.setEnabled(false);
    	
    	final SettingsModelIntegerBounded gapextend=new SettingsModelIntegerBounded(SNPcallNodeModel.CFGKEY_GAPEXTEND, SNPcallNodeModel.DEFAULT_GAPEXTEND, 0, Integer.MAX_VALUE);
    	final SettingsModelDoubleBounded minfrac = new SettingsModelDoubleBounded(SNPcallNodeModel.CFGKEY_MINFRAC,SNPcallNodeModel.DEFAULT_MINFRAC, 0, Double.MAX_VALUE);
    	final SettingsModelIntegerBounded homopoly = new SettingsModelIntegerBounded(SNPcallNodeModel.CFGKEY_HOMOPOLY, SNPcallNodeModel.DEFAULT_HOMOPOLY, 0, Integer.MAX_VALUE);
    	final SettingsModelBoolean noindel = new SettingsModelBoolean(SNPcallNodeModel.CFGKEY_NOINDEL, false);
    	final SettingsModelIntegerBounded skipindel= new SettingsModelIntegerBounded(SNPcallNodeModel.CFGKEY_SKIPINDEL,SNPcallNodeModel.DEFAULT_SKIPINDEL, 0, Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded mingapreads= new SettingsModelIntegerBounded(SNPcallNodeModel.CFGKEY_MINGAPREADS,SNPcallNodeModel.DEFAULT_MINGAPREADS, 0, Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded gapopen = new SettingsModelIntegerBounded(SNPcallNodeModel.CFGKEY_GAPOPEN,SNPcallNodeModel.DEFAULT_GAPOPEN, 0, Integer.MAX_VALUE);
    	
    	createNewGroup("Input Options:");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(
    			SNPcallNodeModel.CFGKEY_ENCODING, false), "Quality is in the Illumina 1.3+ encoding"));
    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(
    			SNPcallNodeModel.CFGKEY_ANAMALOUS, false), "Do not skip anomalous read pairs in variant calling"));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(
    			SNPcallNodeModel.CFGKEY_PROBREALIGN, false), "Disable probabilistic realignment for the computation of base alignment quality (BAQ)"));
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					SNPcallNodeModel.CFGKEY_DOWNGRADE, 
    					SNPcallNodeModel.DEFAULT_DOWNGRADE, 0, Integer.MAX_VALUE),
    					"Coefficient for downgrading mapping quality for reads containing excessive mismatches:", /*step*/ 1));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					SNPcallNodeModel.CFGKEY_MAXREADS, 
    					SNPcallNodeModel.DEFAULT_MAXREADS, 0, Integer.MAX_VALUE),
    					"Maximally reads per input BAM:", /*step*/ 1));
    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(
    			SNPcallNodeModel.CFGKEY_EXTENDBAQ, false), "Extended BAQ computation"));
    	setHorizontalPlacement(false);

    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					SNPcallNodeModel.CFGKEY_CAPMAPQUAL, 
    					SNPcallNodeModel.DEFAULT_CAPMAPQUAL, 0, Integer.MAX_VALUE),
    					"Cap mapping quality at:", /*step*/ 1));
       	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(
    			SNPcallNodeModel.CFGKEY_IGNORERG, false), "Ignore RG tags"));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					SNPcallNodeModel.CFGKEY_MINMAPQUAL, 
    					SNPcallNodeModel.DEFAULT_MINMAPQUAL, 0, Integer.MAX_VALUE),
    					"Minimum mapping quality for an alignment to be used", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					SNPcallNodeModel.CFGKEY_MINBASEQUAL, 
    					SNPcallNodeModel.DEFAULT_MINBASEQUAL, 0, Integer.MAX_VALUE),
    					"Minimum base quality for a base to be considered", /*step*/ 1));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(ifreadfile, "Activate if you want to exclude read groups listed in: "));
    	addDialogComponent(new DialogComponentFileChooser(excludereadfile, "his_id1", 0, false));
    	setHorizontalPlacement(false);
    	
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(ifbedfile, "Activate if you want to restrict to a list of positions (chr pos) or regions (BED):"));
    	addDialogComponent(new DialogComponentFileChooser(bedfile, "his_id", 0, false));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentOptionalString(new SettingsModelOptionalString(SNPcallNodeModel.CFGKEY_PILEUPREGION, "", false), "Region in which pileup is generated:"));
    	
    	createNewGroup("Options for Genotype Likelihood Computation:");
    	addDialogComponent(new DialogComponentNumber(
				minfrac,
				"Minimum fraction of gapped reads for candidates", /*step*/ 0.001));  	
addDialogComponent(new DialogComponentNumber(
				homopoly,
				"Coefficient for modeling homopolymer errors", /*step*/ 1));
addDialogComponent(new DialogComponentBoolean(noindel, "Do not perform INDEL calling"));
addDialogComponent(new DialogComponentNumber(
		gapopen,
		"Phred-scaled gap open sequencing error probability. Reducing x leads to more indel calls:", /*step*/ 1));
addDialogComponent(new DialogComponentNumber(
		gapextend,
	"Phred-scaled gap extension sequencing error probability. Reducing x leads to longer indels", /*step*/ 1));
addDialogComponent(new DialogComponentNumber(
				skipindel,
				"Skip INDEL calling if the average per-sample depth is above:", /*step*/ 1));
addDialogComponent(new DialogComponentNumber(
				mingapreads,
				"Minimum gapped reads for indel candidates:", /*step*/ 1));
addDialogComponent(new DialogComponentOptionalString(new SettingsModelOptionalString(SNPcallNodeModel.CFGKEY_PLATFORMLIST, "", false), "Comma separated list of platforms for indels:"));

    	
    	createNewTab("Filter Options");
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					SNPcallNodeModel.CFGKEY_MINRMS, 
    					SNPcallNodeModel.DEFAULT_MINRMS, 0, Integer.MAX_VALUE),
    					"Minimum RMS mapping quality for SNPs:", /*step*/ 1));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					SNPcallNodeModel.CFGKEY_MINREADDEPTH, 
    					SNPcallNodeModel.DEFAULT_MINREADDEPTH, 0, Integer.MAX_VALUE),
    					"Minimum read depth:", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					SNPcallNodeModel.CFGKEY_MAXREADDEPTH, 
    					SNPcallNodeModel.DEFAULT_MAXREADDEPTH, 0, Integer.MAX_VALUE),
    					"Maximum read depth:", /*step*/ 1));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					SNPcallNodeModel.CFGKEY_MINALTBASE, 
    					SNPcallNodeModel.DEFAULT_MINALTBASE, 0, Integer.MAX_VALUE),
    					"Minimum number of alternate bases:", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					SNPcallNodeModel.CFGKEY_GAPFILTER, 
    					SNPcallNodeModel.DEFAULT_GAPFILTER, 0, Integer.MAX_VALUE),
    					"SNP within x bp around a gap to be filtered:", /*step*/ 1));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					SNPcallNodeModel.CFGKEY_ADJACENTGAPS, 
    					SNPcallNodeModel.DEFAULT_ADJACENTGAPS, 0, Integer.MAX_VALUE),
    					"Window size for filtering adjacent gaps:", /*step*/ 1));
    	createNewGroup("Set p-value thresholds");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelDoubleBounded(
    					SNPcallNodeModel.CFGKEY_STRANDPVAL, 
    					SNPcallNodeModel.DEFAULT_STRANDPVAL, 0, Double.MAX_VALUE),
    					"min P-value for strand bias:", /*step*/ 0.0001));
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					SNPcallNodeModel.CFGKEY_BASEQPVAL, 
    					SNPcallNodeModel.DEFAULT_BASEQPVAL, 0, Integer.MAX_VALUE),
    					"min P-value for baseQ bias [1.0E-]:", /*step*/ 1));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelDoubleBounded(
    					SNPcallNodeModel.CFGKEY_MAPQPVAL, 
    					SNPcallNodeModel.DEFAULT_MAPQPVAL, 0, Double.MAX_VALUE),
    					"min P-value for mapQ bias:", /*step*/ 0.0001));
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelDoubleBounded(
    					SNPcallNodeModel.CFGKEY_ENDDISTPVAL, 
    					SNPcallNodeModel.DEFAULT_ENDDISTPVAL, 0, Double.MAX_VALUE),
    					"min P-value for end distance bias:", /*step*/ 0.0001));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelDoubleBounded(
    					SNPcallNodeModel.CFGKEY_HWEPVAL, 
    					SNPcallNodeModel.DEFAULT_HWEPVAL, 0, Double.MAX_VALUE),
    					"min P-value for HWE (plus F<0):", /*step*/ 0.0001));
    	createNewGroup("Output Options");
    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(
    			SNPcallNodeModel.CFGKEY_PRINTFILTERED, false), "Print filtered variants"));
    	
    	/**Checkboxes for the optional arguments**/
    	
    	ifbedfile.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				bedfile.setEnabled(ifbedfile.getBooleanValue());
				}
		});
    	ifreadfile.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				excludereadfile.setEnabled(ifreadfile.getBooleanValue());
				}
		});

    	noindel.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(noindel.isEnabled()){
					skipindel.setEnabled(!noindel.getBooleanValue());
					gapextend.setEnabled(!noindel.getBooleanValue());
					mingapreads.setEnabled(!noindel.getBooleanValue());
					gapopen.setEnabled(!noindel.getBooleanValue());
				}


				}
		});
                    
    }
}

