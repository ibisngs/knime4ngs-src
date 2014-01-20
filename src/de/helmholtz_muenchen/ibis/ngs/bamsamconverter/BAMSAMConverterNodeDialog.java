package de.helmholtz_muenchen.ibis.ngs.bamsamconverter;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "BAMSAMConverter" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class BAMSAMConverterNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring the BAMSAMConverter node.
     */
	
	
	final SettingsModelString path2samtools = new SettingsModelString(
			BAMSAMConverterNodeModel.CFGKEY_PATH2SAMTOOLS,"");
	final SettingsModelString method = new SettingsModelString(
			BAMSAMConverterNodeModel.CFGKEY_METHOD,"");
	final SettingsModelString infile = new SettingsModelString(
			BAMSAMConverterNodeModel.CFGKEY_INFILE,"");
//	final SettingsModelBoolean outputbam = new SettingsModelBoolean(
//			BAMSAMConverterNodeModel.CFGKEY_outputbam, true);
	final SettingsModelBoolean printsamheader = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_printsamheader, true);
	final SettingsModelBoolean printsamheaderonly = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_printsamheaderonly, false);
	final SettingsModelBoolean uncompressedbam = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_uncompressedbam, false);
	final SettingsModelBoolean fastcompression = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_fastcompression, false);
	final SettingsModelBoolean outhexflag = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_outhexflag, false);
	final SettingsModelBoolean outstringflag = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_outstringflag, false);
	final SettingsModelBoolean printmatchingrecordsonly = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_printmatchingrecordsonly, false);
	final SettingsModelString bedfile = new SettingsModelString(
			BAMSAMConverterNodeModel.CFGKEY_bedfile, "");
	final SettingsModelString refnamelist = new SettingsModelString(
			BAMSAMConverterNodeModel.CFGKEY_refnamelist, "");
	final SettingsModelString refseqfile = new SettingsModelString(
			BAMSAMConverterNodeModel.CFGKEY_refseqfile, "");
	final SettingsModelString listofreads = new SettingsModelString(
			BAMSAMConverterNodeModel.CFGKEY_listofreads, "");	
	final SettingsModelIntegerBounded requiredflag = new SettingsModelIntegerBounded(
			BAMSAMConverterNodeModel.CFGKEY_requiredflag,0,0,Integer.MAX_VALUE);
	final SettingsModelIntegerBounded filteringflag = new SettingsModelIntegerBounded(
			BAMSAMConverterNodeModel.CFGKEY_filteringflag,0,0,Integer.MAX_VALUE);
	final SettingsModelIntegerBounded minmapqual = new SettingsModelIntegerBounded(
			BAMSAMConverterNodeModel.CFGKEY_minmapqual,0,0,Integer.MAX_VALUE);
	final SettingsModelString onlylibraryreads = new SettingsModelString(
			BAMSAMConverterNodeModel.CFGKEY_onlylibraryreads, "");
	final SettingsModelString onlygroupreads = new SettingsModelString(
			BAMSAMConverterNodeModel.CFGKEY_onlygroupreads, "");
	
	//Checkboxes
	final SettingsModelBoolean usebedfile = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_USEBEDFILE, false);
	final SettingsModelBoolean userefnamelist = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_USEREFNAMELIST, false);
	final SettingsModelBoolean userefseqfile = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_USEREFSEQFILE, false);
	final SettingsModelBoolean userequiredflag = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_USEREQUIREDFLAG, false);
	final SettingsModelBoolean usefilteringflag = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_USEFILTERINGFLAG, false);
	final SettingsModelBoolean useonlylibraryreads = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_USEONLYLIBRARYREADS, false);
	final SettingsModelBoolean useonlygroupreads = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_USEONLYGROUPREADS, false);
	final SettingsModelBoolean uselistofreads = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_USELISTOFREADS, false);
	 final SettingsModelBoolean checkbamindex = new SettingsModelBoolean(
			BAMSAMConverterNodeModel.CFGKEY_CHECKBAMINDEX, true);
/**
 * 	/**Usage:   samtools view [options] <in.bam>|<in.sam> [region1 [...]]

	Options: -b       output BAM
	         -h       print header for the SAM output
	         -H       print header only (no alignments)
	         -S       input is SAM
	         -u       uncompressed BAM output (force -b)
	         -1       fast compression (force -b)
	         -x       output FLAG in HEX (samtools-C specific)
	         -X       output FLAG in string (samtools-C specific)
	         -c       print only the count of matching records
	         -L FILE  output alignments overlapping the input BED FILE [null]
	         -t FILE  list of reference names and lengths (force -S) [null]
	         -T FILE  reference sequence file (force -S) [null]
	         -o FILE  output file name [stdout]
	         -R FILE  list of read groups to be outputted [null]


	 * 
	 */
    protected BAMSAMConverterNodeDialog() {
    	
 //   	outputbam.setBooleanValue(true);
    	printsamheader.setBooleanValue(true);
    	
    	
    	
    	createNewGroup("Path to Samtools");
    	addDialogComponent(new DialogComponentFileChooser(path2samtools, "his_bamsam_id5", 0, ""));
    	createNewGroup("");
    	addDialogComponent(new DialogComponentStringSelection(method,"Select utility","SAM->BAM","BAM->SAM"));
    	createNewGroup("Select file to convert");
    	addDialogComponent(new DialogComponentFileChooser(infile,"his_bamsam_ID0",".sam",".bam"));
	
    	createNewGroup("Options");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(checkbamindex,"Index BAM file"));	
    	//addDialogComponent(new DialogComponentBoolean(outputbam,"Output in BAM format"));
    	addDialogComponent(new DialogComponentBoolean(uncompressedbam,"Uncompressed BAM output"));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(fastcompression,"Fast compression"));
    	addDialogComponent(new DialogComponentBoolean(printsamheader,"Print header for the SAM output"));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(printsamheaderonly,"Print header only (no alignments)"));
    	addDialogComponent(new DialogComponentBoolean(printmatchingrecordsonly,"Print only the count of matching records"));
    	setHorizontalPlacement(false);

    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(outhexflag,"Output FLAG in HEX"));
    	addDialogComponent(new DialogComponentBoolean(outstringflag,"Output FLAG in string"));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usebedfile,"Output alignments overlapping the input BED file"));
    	addDialogComponent(new DialogComponentFileChooser(bedfile, "his_bamsam_id1", 0, ""));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(userefnamelist,"List of reference names and lengths"));
    	addDialogComponent(new DialogComponentFileChooser(refnamelist, "his_bamsam_id2", 0, ""));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(userefseqfile,"Reference sequence file"));
    	addDialogComponent(new DialogComponentFileChooser(refseqfile, "his_bamsam_id3", 0, ""));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(uselistofreads,"List of read groups to be outputted"));
    	addDialogComponent(new DialogComponentFileChooser(listofreads, "his_bamsam_id4", 0, ""));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(userequiredflag,"Required flag"));
    	addDialogComponent(new DialogComponentNumber(requiredflag, "",1));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usefilteringflag,"Filtering flag"));
    	addDialogComponent(new DialogComponentNumber(filteringflag, "",1));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentNumber(minmapqual, "Minimum mapping quality",1));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(useonlylibraryreads,"Only output reads in library:"));
    	addDialogComponent(new DialogComponentString(onlylibraryreads,""));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(useonlygroupreads,"Only output reads in read group:"));
    	addDialogComponent(new DialogComponentString(onlygroupreads, ""));
    	setHorizontalPlacement(false);
    	
    	
    	usebedfile.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				bedfile.setEnabled(usebedfile.getBooleanValue());
				}
		});
    	userefnamelist.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				refnamelist.setEnabled(userefnamelist.getBooleanValue());
				}
		});
    	userefseqfile.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				refseqfile.setEnabled(userefseqfile.getBooleanValue());
				}
		});
    	uselistofreads.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				listofreads.setEnabled(uselistofreads.getBooleanValue());
				}
		});
    	
    	
    	userequiredflag.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				requiredflag.setEnabled(userequiredflag.getBooleanValue());
				}
		});
    	usefilteringflag.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				filteringflag.setEnabled(usefilteringflag.getBooleanValue());
				}
		});
    	useonlygroupreads.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				onlygroupreads.setEnabled(useonlygroupreads.getBooleanValue());
				}
		});
    	useonlylibraryreads.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				onlylibraryreads.setEnabled(useonlylibraryreads.getBooleanValue());
				}
		});
    	
    	method.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(method.getStringValue().equals("SAM->BAM")){
//			    	outputbam.setBooleanValue(true);
//			    	outputbam.setEnabled(true);
			    	printsamheader.setBooleanValue(false);
			    	printsamheader.setEnabled(false);
			    	printsamheaderonly.setBooleanValue(false);
			    	printsamheaderonly.setEnabled(false);
			    	fastcompression.setEnabled(true);
			    	uncompressedbam.setEnabled(true);
			    	userefnamelist.setEnabled(true);
			    	userefseqfile.setEnabled(true);
			    	checkbamindex.setEnabled(true);

				}else{
//			    	outputbam.setBooleanValue(false);
//			    	outputbam.setEnabled(false);
			    	printsamheaderonly.setEnabled(true);
			    	printsamheader.setEnabled(true);
			    	fastcompression.setEnabled(false);
			    	uncompressedbam.setEnabled(false);
			    	userefseqfile.setEnabled(false);
			    	userefnamelist.setEnabled(false);
			    	checkbamindex.setEnabled(false);
				}
				}
		});
    	
    	uncompressedbam.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(uncompressedbam.getBooleanValue()){
//					outputbam.setBooleanValue(true);
				}
				}
		});
    	fastcompression.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(fastcompression.getBooleanValue()){
//					outputbam.setBooleanValue(true);
				}
				}
		});
    	
    }
}

