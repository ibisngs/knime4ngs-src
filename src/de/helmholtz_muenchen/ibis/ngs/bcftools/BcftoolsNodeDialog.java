package de.helmholtz_muenchen.ibis.ngs.bcftools;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "Bcftools" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Max
 */
public class BcftoolsNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring the Bcftools node.
     */
    protected BcftoolsNodeDialog() {

    	final SettingsModelString bcfmethod = new SettingsModelString(BcftoolsNodeModel.CFGKEY_BCFMETHOD,"");
    	final SettingsModelString catinfile = new SettingsModelString(BcftoolsNodeModel.CFGKEY_CATINFILE,"");
    	final SettingsModelString infile = new SettingsModelString(BcftoolsNodeModel.CFGKEY_INFILE,"");
    	final SettingsModelString vcfsampleheader = new SettingsModelString(BcftoolsNodeModel.CFGKEY_VCFSAMPLEHEADER, "");
    	catinfile.setEnabled(false);
    	vcfsampleheader.setEnabled(false);
    	
    	/**
    	 *Models for Bcftools call
    	 */
    	final SettingsModelBoolean outbcf = new SettingsModelBoolean(BcftoolsNodeModel.CFGKEY_OUTBCF, false);
    	final SettingsModelBoolean outuncompressedbcf = new SettingsModelBoolean(BcftoolsNodeModel.CFGKEY_OUTUNCOMPRESSEDBCF, false);
    	final SettingsModelBoolean ifbedfile =  new SettingsModelBoolean(BcftoolsNodeModel.CFGKEY_IFBEDFILE, false);
    	final SettingsModelBoolean ifseqdic =  new SettingsModelBoolean(BcftoolsNodeModel.CFGKEY_IFSEQDIC, false);
    	final SettingsModelBoolean ifsamplelist =  new SettingsModelBoolean(BcftoolsNodeModel.CFGKEY_IFSAMPLELIST, false);
    	final SettingsModelString bedfile = new SettingsModelString(BcftoolsNodeModel.CFGKEY_BEDFILE, "");
    	final SettingsModelString seqdic = new SettingsModelString(BcftoolsNodeModel.CFGKEY_SEQDIC, "");
    	final SettingsModelString samplelist = new SettingsModelString(BcftoolsNodeModel.CFGKEY_SAMPLELIST, "");
    	bedfile.setEnabled(false);
    	seqdic.setEnabled(false);
    	samplelist.setEnabled(false);
    	
    	
    	final SettingsModelBoolean snpcalling = new SettingsModelBoolean(BcftoolsNodeModel.CFGKEY_SNPCALLING, false);
    	final SettingsModelBoolean callgenotype = new SettingsModelBoolean(BcftoolsNodeModel.CFGKEY_CALLGENOTYPE, false);
    	final SettingsModelDoubleBounded skiploci = new SettingsModelDoubleBounded(BcftoolsNodeModel.CFGKEY_SAMPLECOVERAGE, BcftoolsNodeModel.DEFAULT_SAMPLECOVERAGE, 0, Double.MAX_VALUE);
    	final SettingsModelBoolean outvariantsonly = new SettingsModelBoolean(BcftoolsNodeModel.CFGKEY_OUTVARIANTSONLY, false);
    	final SettingsModelBoolean likelihoodana = new SettingsModelBoolean(BcftoolsNodeModel.CFGKEY_LIKELIHOODANA, false);
    	final SettingsModelBoolean skipindel = new SettingsModelBoolean(BcftoolsNodeModel.CFGKEY_SKIPINDEL, false);
    	final SettingsModelDoubleBounded indelsubratio = new SettingsModelDoubleBounded(BcftoolsNodeModel.CFGKEY_INDELSUBRATIO, BcftoolsNodeModel.DEFAULT_INDELSUBRATIO, -Double.MAX_VALUE, Double.MAX_VALUE);
    	final SettingsModelString constrainedcalling = new SettingsModelString(BcftoolsNodeModel.CFGKEY_CONSTRAINEDCALLING,"No constrains");
    	skiploci.setEnabled(false);
    	
    	//Main Tab 
    	createNewGroup("Path to Bcftools");
//    	addDialogComponent(new DialogComponentFileChooser(use_FILE, "Testing", 0, true));
    	addDialogComponent(new DialogComponentFileChooser(new SettingsModelString(BcftoolsNodeModel.CFGKEY_PATH2BCFTOOLS, ""),"his_bcft_ID5",""));
    	createNewGroup("");
    	addDialogComponent(new DialogComponentStringSelection(bcfmethod,"Select method","call","index", "concat","reheader"));
    	createNewGroup("Select Inputfile for selected method");
    	addDialogComponent(new DialogComponentFileChooser(infile,"his_bcft_ID3",".bcf",".vcf"));
    	createNewGroup("Select second BCF file (concat only)");
    	addDialogComponent(new DialogComponentFileChooser(catinfile,"his_bcft_ID4",".bcf"));
    	createNewGroup("Select file that includes the new sample names");
    	addDialogComponent(new DialogComponentFileChooser(vcfsampleheader,"hisID_bcft_vcfsampleheader",""));
    	
    	
    	
    	createNewTab("View In/Out Options");
    	createNewGroup("Input/Output Options");
//    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(
//    			BcftoolsNodeModel.CFGKEY_INISVCF, false), "Input is VCF"));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(ifseqdic, "Sequence dictionary for VCF->BCF conversion:"));
    	addDialogComponent(new DialogComponentFileChooser(seqdic, "his_bcft_ID",""));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(ifbedfile, "List of sites (chr pos) or regions (BED) to output:"));
    	addDialogComponent(new DialogComponentFileChooser(bedfile, "his_bcft_ID1",""));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(ifsamplelist, "List of samples to use:"));
    	addDialogComponent(new DialogComponentFileChooser(samplelist, "his_bcft_ID2",""));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(
    			BcftoolsNodeModel.CFGKEY_KEEPALLELES, false), "Keep all possible alternate alleles at variant sites"));
    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(
    			BcftoolsNodeModel.CFGKEY_CALCLD, false), "Calculate LD for adjacent sites"));
    	setHorizontalPlacement(false);
    	
    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(
    			BcftoolsNodeModel.CFGKEY_PLGENERATE, false), "PL generated by r921 or before (which generate old ordering)"));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(
    			BcftoolsNodeModel.CFGKEY_SURPESSGENOTYPE, false), "Suppress all individual genotype information"));
    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(
    			BcftoolsNodeModel.CFGKEY_SKIPREF, false), "Skip sites where REF is not A/C/G/T"));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(outbcf, "Output BCF instead of VCF"));
    	addDialogComponent(new DialogComponentBoolean(outuncompressedbcf, "Uncompressed BCF output"));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(
    			BcftoolsNodeModel.CFGKEY_OUTQCALL, false), "Output the QCALL likelihood format"));


    	createNewTab("Calling Options");
    	createNewGroup("Consensus/variant calling options");
    	//createNewGroup("Consensus/variant calling options:");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(snpcalling, "Call SNPs"));
    	addDialogComponent(new DialogComponentBoolean(likelihoodana, "Likelihood based analyses"));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(callgenotype, "Call genotypes at variant sites"));
    	addDialogComponent(new DialogComponentBoolean(outvariantsonly, "Output potential variant sites only"));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentStringSelection(
    			new SettingsModelString(BcftoolsNodeModel.CFGKEY_TYPEOFPRIOR,BcftoolsNodeModel.DEFAULT_TYPEOFPRIOR),
    			"Type of prior","full", "cond2", "flat"));
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelDoubleBounded(
    					BcftoolsNodeModel.CFGKEY_VARIANTIF, 
    					BcftoolsNodeModel.DEFAULT_VARIANTIF, 0, Double.MAX_VALUE),
    					"Variant if P(ref|D)<:", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(
    			skiploci,
    					"Skip loci if fraction of samples covered is less than:", /*step*/ 1));
    	addDialogComponent(new DialogComponentStringSelection(
    			constrainedcalling,
    			"Constrained calling","No constrains","pair", "trioauto", "trioxd","trioxs"));
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelDoubleBounded(
    					BcftoolsNodeModel.CFGKEY_MUTATIONRATE, 
    					BcftoolsNodeModel.DEFAULT_MUTATIONRATE, 0, Double.MAX_VALUE),
    					"Scaled substitution mutation rate", /*step*/ 1));
    	addDialogComponent(new DialogComponentBoolean(skipindel, "Skip indels"));
    	addDialogComponent(new DialogComponentNumber(indelsubratio,"Indel-to-substitution ratio", /*step*/ 1));

    	
    	createNewGroup("Contrast calling and association test options:");
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					BcftoolsNodeModel.CFGKEY_NUMBEROFGRPSAMPLES, 
    					BcftoolsNodeModel.DEFAULT_NUMBEROFGRPSAMPLES, 0, Integer.MAX_VALUE),
    					"Number of group-1 samples", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelDoubleBounded(
    					BcftoolsNodeModel.CFGKEY_POSTERIORICON, 
    					BcftoolsNodeModel.DEFAULT_POSTERIORICON, 0, Double.MAX_VALUE),
    					"Posterior contrast for P(ref|D)<0.5 and LRT<", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					BcftoolsNodeModel.CFGKEY_NUMBEROFPERMUTAIONS, 
    					BcftoolsNodeModel.DEFAULT_NUMBEROFPERMUTAIONS, 0, Integer.MAX_VALUE),
    					"Number of permutations for association testing", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelDoubleBounded(
    					BcftoolsNodeModel.CFGKEY_ONLYPERMUTAIONS, 
    					BcftoolsNodeModel.DEFAULT_ONLYPERMUTAIONS, 0, Double.MAX_VALUE),
    					"Only perform permutations for P(chi^2)<", /*step*/ 1));
    	
    	
    	
    	//Main method chooser
    	bcfmethod.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(bcfmethod.getStringValue().equals("view")){
					setEnabled(true,"View In/Out Options");
			    	setEnabled(true,"Call Options");
				}else{
					setEnabled(false,"View In/Out Options");
			    	setEnabled(false,"Calling Options");
				}
				if(bcfmethod.getStringValue().equals("concat")){
					catinfile.setEnabled(true);
				}else{
					catinfile.setEnabled(false);
				}
				if(bcfmethod.getStringValue().equals("reheader")){
					vcfsampleheader.setEnabled(true);
				}else{
					vcfsampleheader.setEnabled(false);
				}
				
				}
		});
    	
    	
    	
    	outuncompressedbcf.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(outuncompressedbcf.getBooleanValue()){
					outbcf.setBooleanValue(true);
				}
				}
		});
    	outbcf.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(!outbcf.getBooleanValue()){
					outuncompressedbcf.setBooleanValue(false);
				}
				}
		});
    	ifbedfile.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				bedfile.setEnabled(ifbedfile.getBooleanValue());
				}
		});
    	ifseqdic.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				seqdic.setEnabled(ifseqdic.getBooleanValue());
				}
		});
    	ifsamplelist.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				samplelist.setEnabled(ifsamplelist.getBooleanValue());
				}
		});
    	snpcalling.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(snpcalling.getBooleanValue()){
					likelihoodana.setBooleanValue(true);
				}else{
					callgenotype.setBooleanValue(false);
					outvariantsonly.setBooleanValue(false);
				}
				}
		});
    	likelihoodana.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(!likelihoodana.getBooleanValue()){
					snpcalling.setBooleanValue(false);
				}
				}
		});
    	callgenotype.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(callgenotype.getBooleanValue()){
					snpcalling.setBooleanValue(true);
				}
				}
		});
    	outvariantsonly.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(outvariantsonly.getBooleanValue()){
					snpcalling.setBooleanValue(true);
				}
				skiploci.setEnabled(outvariantsonly.getBooleanValue());
				}
		});
    	skipindel.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				indelsubratio.setEnabled(!skipindel.getBooleanValue());
				}
		});
    	constrainedcalling.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(!constrainedcalling.getStringValue().equals("pair")&&!constrainedcalling.getStringValue().equals("No constrains")){
					ifsamplelist.setBooleanValue(true);
				}
				}
		}); 	
 
    	
    }
}

