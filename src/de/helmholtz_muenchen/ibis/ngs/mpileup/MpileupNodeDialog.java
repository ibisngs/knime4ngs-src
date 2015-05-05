package de.helmholtz_muenchen.ibis.ngs.mpileup;

import javax.swing.JFileChooser;
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

import de.helmholtz_muenchen.ibis.ngs.bwa.BWANodeModel;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeModel;

/**
 * <code>NodeDialog</code> for the "Mpileup" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Max
 */
public class MpileupNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring the Mpileup node.
     */
    protected MpileupNodeDialog() {
    	/*Mpileup
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

    Output options:

           -D           output per-sample DP in BCF (require -g/-u)
           -g           generate BCF output (genotype likelihoods)
           -O           output base positions on reads (disabled by -g/-u)
           -s           output mapping quality (disabled by -g/-u)
           -S           output per-sample strand bias P-value in BCF (require -g/-u)
           -u           generate uncompress BCF output

    ##SNP Call specific parameters
           -e INT       Phred-scaled gap extension seq error probability [20]
           -F FLOAT     minimum fraction of gapped reads for candidates [0.002]
           -h INT       coefficient for homopolymer errors [100]
           -I           do not perform indel calling
           -L INT       max per-sample depth for INDEL calling [250]
           -m INT       minimum gapped reads for indel candidates [1]
           -o INT       Phred-scaled gap open sequencing error probability [40]
           -P STR       comma separated list of platforms for indels [all]

    	 */
        final SettingsModelString SAM = new SettingsModelString(MpileupNodeModel.CFGKEY_SAM_PATH, "");
    	final SettingsModelString bedfile=new SettingsModelString(MpileupNodeModel.CFGKEY_BEDFILE,null);
		final SettingsModelBoolean ifbedfile = new SettingsModelBoolean(MpileupNodeModel.CFGKEY_IFBEDFILE, false);
    	bedfile.setEnabled(false);
    	final SettingsModelString excludereadfile=new SettingsModelString(MpileupNodeModel.CFGKEY_EXCLUDEREADS,null);
		final SettingsModelBoolean ifreadfile = new SettingsModelBoolean(MpileupNodeModel.CFGKEY_IFREADFILE, false);
    	excludereadfile.setEnabled(false);
    	
    	final SettingsModelBoolean outpersample = new SettingsModelBoolean(MpileupNodeModel.CFGKEY_OUTPERSAMPLE, false);
    	final SettingsModelBoolean  bcfoutput = 	new SettingsModelBoolean(MpileupNodeModel.CFGKEY_BCFOUTPUT, true);
    	final SettingsModelBoolean outbasepositions = new SettingsModelBoolean(MpileupNodeModel.CFGKEY_OUTBASEPOSITIONS, false);
    	outbasepositions.setEnabled(false);
    	final SettingsModelBoolean outmapqual = new SettingsModelBoolean(MpileupNodeModel.CFGKEY_OUTMAPQUAL, false);
    	outmapqual.setEnabled(false);
    	final SettingsModelBoolean strandbiaspval = new SettingsModelBoolean(MpileupNodeModel.CFGKEY_STRANDBIASPVAL, false);
    	final SettingsModelBoolean uncompressedbcf = new SettingsModelBoolean(MpileupNodeModel.CFGKEY_UNCOMPRESSEDBCF, false);
    	
    	final SettingsModelIntegerBounded gapextend=new SettingsModelIntegerBounded(MpileupNodeModel.CFGKEY_GAPEXTEND, MpileupNodeModel.DEFAULT_GAPEXTEND, 0, Integer.MAX_VALUE);
    	final SettingsModelDoubleBounded minfrac = new SettingsModelDoubleBounded(MpileupNodeModel.CFGKEY_MINFRAC,MpileupNodeModel.DEFAULT_MINFRAC, 0, Double.MAX_VALUE);
    	final SettingsModelIntegerBounded homopoly = new SettingsModelIntegerBounded(MpileupNodeModel.CFGKEY_HOMOPOLY, MpileupNodeModel.DEFAULT_HOMOPOLY, 0, Integer.MAX_VALUE);
    	final SettingsModelBoolean noindel = new SettingsModelBoolean(MpileupNodeModel.CFGKEY_NOINDEL, false);
    	final SettingsModelIntegerBounded skipindel= new SettingsModelIntegerBounded(MpileupNodeModel.CFGKEY_SKIPINDEL,MpileupNodeModel.DEFAULT_SKIPINDEL, 0, Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded mingapreads= new SettingsModelIntegerBounded(MpileupNodeModel.CFGKEY_MINGAPREADS,MpileupNodeModel.DEFAULT_MINGAPREADS, 0, Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded gapopen = new SettingsModelIntegerBounded(MpileupNodeModel.CFGKEY_GAPOPEN,MpileupNodeModel.DEFAULT_GAPOPEN, 0, Integer.MAX_VALUE);
    	
    	final SettingsModelBoolean usefaidxfile = new SettingsModelBoolean(
    			MpileupNodeModel.CFGKEY_USEFAIDX,true);
    	final SettingsModelBoolean ifindex = new SettingsModelBoolean(
    			MpileupNodeModel.CFGKEY_IFINDEX,true);
    	
    	final SettingsModelOptionalString platformlist = new SettingsModelOptionalString(MpileupNodeModel.CFGKEY_PLATFORMLIST, "", false);
    	
    	createNewGroup("Path to samtools software");
    	addDialogComponent(new DialogComponentFileChooser(SAM, "sam_mpileup", 0, ""));

    	createNewGroup("Input Options:");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usefaidxfile, "Use indexed reference genome"));
    	addDialogComponent(new DialogComponentBoolean(ifindex, "Index reference genome first"));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(
    			MpileupNodeModel.CFGKEY_ENCODING, false), "Quality is in the Illumina 1.3+ encoding"));
    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(
    			MpileupNodeModel.CFGKEY_ANAMALOUS, false), "Do not skip anomalous read pairs in variant calling"));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(
    			MpileupNodeModel.CFGKEY_PROBREALIGN, false), "Disable probabilistic realignment for the computation of base alignment quality (BAQ)"));
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					MpileupNodeModel.CFGKEY_DOWNGRADE, 
    					MpileupNodeModel.DEFAULT_DOWNGRADE, 0, Integer.MAX_VALUE),
    					"Coefficient for downgrading mapping quality for reads containing excessive mismatches:", /*step*/ 1));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					MpileupNodeModel.CFGKEY_MAXREADS, 
    					MpileupNodeModel.DEFAULT_MAXREADS, 0, Integer.MAX_VALUE),
    					"Maximally reads per input BAM:", /*step*/ 1));
    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(
    			MpileupNodeModel.CFGKEY_EXTENDBAQ, false), "Extended BAQ computation"));
    	setHorizontalPlacement(false);

    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					MpileupNodeModel.CFGKEY_CAPMAPQUAL, 
    					MpileupNodeModel.DEFAULT_CAPMAPQUAL, 0, Integer.MAX_VALUE),
    					"Cap mapping quality at:", /*step*/ 1));
       	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(
    			MpileupNodeModel.CFGKEY_IGNORERG, false), "Ignore RG tags"));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					MpileupNodeModel.CFGKEY_MINMAPQUAL, 
    					MpileupNodeModel.DEFAULT_MINMAPQUAL, 0, Integer.MAX_VALUE),
    					"Minimum mapping quality for an alignment to be used", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					MpileupNodeModel.CFGKEY_MINBASEQUAL, 
    					MpileupNodeModel.DEFAULT_MINBASEQUAL, 0, Integer.MAX_VALUE),
    					"", /*step*/ 1));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(ifreadfile, "Activate if you want to exclude read groups listed in: "));
    	addDialogComponent(new DialogComponentFileChooser(excludereadfile, "his_id1", 0, false));
    	setHorizontalPlacement(false);
    	
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(ifbedfile, "Activate if you want to restrict to a list of positions (chr pos) or regions (BED):"));
    	addDialogComponent(new DialogComponentFileChooser(bedfile, "his_id", 0, false));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentOptionalString(new SettingsModelOptionalString(MpileupNodeModel.CFGKEY_PILEUPREGION, "", false), "Region in which pileup is generated:"));
    	
    	
    	createNewGroup("Output Options");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(bcfoutput, "Generate BCF output (genotype likelihoods)"));
    	addDialogComponent(new DialogComponentBoolean(uncompressedbcf, "Generate uncompress BCF output"));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(outpersample, "Output per-sample DP in BCF"));
    	addDialogComponent(new DialogComponentBoolean(strandbiaspval, "Output per-sample strand bias P-value in BCF"));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(outbasepositions, "Output base positions on reads"));
    	addDialogComponent(new DialogComponentBoolean(outmapqual, "Output mapping quality"));
    	setHorizontalPlacement(false);

    	
    	
    	createNewTab("Options for Genotype Likelihood Computation:"); 
    	addDialogComponent(new DialogComponentNumber(
    					minfrac,
    					"Minimum fraction of gapped reads for candidates", /*step*/ 0.001));  	
    	addDialogComponent(new DialogComponentNumber(
    					homopoly,
    					"Coefficient for modeling homopolymer errors", /*step*/ 1));
    	addDialogComponent(new DialogComponentBoolean(noindel, "Do not perform INDEL calling"));
    	addDialogComponent(new DialogComponentNumber(
				gapextend,
			"Phred-scaled gap extension sequencing error probability", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(
				gapopen,
				"Phred-scaled gap open sequencing error probability", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(
    					skipindel,
    					"Skip INDEL calling if the average per-sample depth is above:", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(
    					mingapreads,
    					"Minimum gapped reads for indel candidates:", /*step*/ 1));
    	addDialogComponent(new DialogComponentOptionalString(platformlist, "Comma separated list of platforms for indels:"));
    
    	
    	
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
       	
    	uncompressedbcf.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(bcfoutput.getBooleanValue() || uncompressedbcf.getBooleanValue()){
					outpersample.setEnabled(true);
					strandbiaspval.setEnabled(true);
					outbasepositions.setEnabled(false);
					outmapqual.setEnabled(false);
					gapextend.setEnabled(true);
					minfrac.setEnabled(true);
					homopoly.setEnabled(true);
					noindel.setEnabled(true);
					skipindel.setEnabled(true);
					mingapreads.setEnabled(true);
					gapopen.setEnabled(true);
					
				}
				if(!uncompressedbcf.getBooleanValue()&&!bcfoutput.getBooleanValue()){
					outbasepositions.setEnabled(true);
					outmapqual.setEnabled(true);
					outpersample.setEnabled(false);
					strandbiaspval.setEnabled(false);
					gapextend.setEnabled(false);
					minfrac.setEnabled(false);
					homopoly.setEnabled(false);
					noindel.setEnabled(false);
					skipindel.setEnabled(false);
					mingapreads.setEnabled(false);
					gapopen.setEnabled(false);
					
				}			
				}
		});
    	bcfoutput.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(bcfoutput.getBooleanValue() || uncompressedbcf.getBooleanValue()){
					outpersample.setEnabled(true);
					strandbiaspval.setEnabled(true);
					outbasepositions.setEnabled(false);
					outmapqual.setEnabled(false);
					gapextend.setEnabled(true);
					minfrac.setEnabled(true);
					homopoly.setEnabled(true);
					noindel.setEnabled(true);
					skipindel.setEnabled(true);
					mingapreads.setEnabled(true);
					gapopen.setEnabled(true);
				}
				if(!uncompressedbcf.getBooleanValue()&&!bcfoutput.getBooleanValue()){
					outbasepositions.setEnabled(true);
					outmapqual.setEnabled(true);
					outpersample.setEnabled(false);
					strandbiaspval.setEnabled(false);
					gapextend.setEnabled(false);
					minfrac.setEnabled(false);
					homopoly.setEnabled(false);
					noindel.setEnabled(false);
					skipindel.setEnabled(false);
					mingapreads.setEnabled(false);
					gapopen.setEnabled(false);
				}			
				}
		});
    	
    	noindel.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(noindel.isEnabled()){
					skipindel.setEnabled(!noindel.getBooleanValue());
					gapextend.setEnabled(!noindel.getBooleanValue());
					mingapreads.setEnabled(!noindel.getBooleanValue());
					gapopen.setEnabled(!noindel.getBooleanValue());
					platformlist.setEnabled(!noindel.getBooleanValue());
					if(noindel.getBooleanValue()){
						platformlist.setIsActive(false);
					}
				}


				}
		});
    	
    	bcfoutput.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
					if(bcfoutput.getBooleanValue()){
						uncompressedbcf.setEnabled(false);
					}else{
						uncompressedbcf.setEnabled(true);
					}
				}
		});
    	uncompressedbcf.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
					if(uncompressedbcf.getBooleanValue()){
						bcfoutput.setEnabled(false);
					}else{
						bcfoutput.setEnabled(true);
					}
				}
		});
    	ifindex.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
					if(ifindex.getBooleanValue()){
						usefaidxfile.setBooleanValue(true);
					}
				}
		});
    	usefaidxfile.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
					if(!usefaidxfile.getBooleanValue()){
						ifindex.setBooleanValue(false);
					}
				}
		});
    	
    	
    }
}

