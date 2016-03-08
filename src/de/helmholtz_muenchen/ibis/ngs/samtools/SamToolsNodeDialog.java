package de.helmholtz_muenchen.ibis.ngs.samtools;

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
 */
public class SamToolsNodeDialog extends HTExecutorNodeDialog {

	final SettingsModelString utility = new SettingsModelString(SamToolsNodeModel.CFGKEY_UTILITY, "");

	final SettingsModelString samtools = new SettingsModelString(SamToolsNodeModel.CFGKEY_SAMTOOLS, "");
	final SettingsModelString bamfile = new SettingsModelString(SamToolsNodeModel.CFGKEY_BAMFILE, "");
	final SettingsModelString refseqfile = new SettingsModelString(SamToolsNodeModel.CFGKEY_REFSEQFILE, "");
	
	//calmd
	final SettingsModelBoolean changeIdentBases = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_CHANGEIDENTBASES, false);
	final SettingsModelBoolean useCompression = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_USECOMPRESSION, false);	
	final SettingsModelString compression = new SettingsModelString(SamToolsNodeModel.CFGKEY_COMPRESSION, "");
	final SettingsModelBoolean inputIsSAM = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_INPUTISSAM, false);
	final SettingsModelBoolean modifyQual = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MODIFYQUAL, false);
	final SettingsModelBoolean bqTag = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_BQTAG, false);
	final SettingsModelBoolean extendedBAQ = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_EXTENDEDBAQ, false);
	//final SettingsModelBoolean doCapMapQual = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_DOCAPMAPQUAL, false);
	//final SettingsModelIntegerBounded capMapQual = new SettingsModelIntegerBounded(SamToolsNodeModel.CFGKEY_CAPMAPQUAL, 1, 1, Integer.MAX_VALUE);
	//rmdup
	final SettingsModelBoolean removeDup = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_REMOVEDUP, false);
	final SettingsModelBoolean treatPE = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_TREATPE, false);
	//cat
	final SettingsModelBoolean useHeaderSAM = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_USEHEADERSAM, false);
	final SettingsModelString headerSAM = new SettingsModelString(SamToolsNodeModel.CFGKEY_HEADERSAM, "");
	final SettingsModelString inBAM1 = new SettingsModelString(SamToolsNodeModel.CFGKEY_INBAM1, "");
	//final SettingsModelString inBAM2 = new SettingsModelString(SamToolsNodeModel.CFGKEY_INBAM2, "");
	//reheader
	//final SettingsModelString rehInBAM = new SettingsModelString(SamToolsNodeModel.CFGKEY_REHINBAM, "");
	final SettingsModelString rehInSAM = new SettingsModelString(SamToolsNodeModel.CFGKEY_REHINSAM, "");
	//merge
	SettingsModelString minbam1 = new SettingsModelString(SamToolsNodeModel.CFGKEY_MINBAM1, "");
	//SettingsModelString minbam2 = new SettingsModelString(SamToolsNodeModel.CFGKEY_MINBAM2, "");
	final SettingsModelBoolean mcompression = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MCOMPRESSION, false);
	final SettingsModelBoolean mforce = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MFORCE, false);
	final SettingsModelBoolean usemhfile = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_USEMHFILE, false);
	final SettingsModelString mhfile = new SettingsModelString(SamToolsNodeModel.CFGKEY_MHFILE, "");
	final SettingsModelBoolean msorted = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MSORTED, false);
	final SettingsModelString mregion = new SettingsModelString(SamToolsNodeModel.CFGKEY_MREGION, "");
	final SettingsModelBoolean usemregion = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_USEMREGION, false);
	final SettingsModelBoolean mrgtag = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MRGTAG, false);
	final SettingsModelBoolean muncompressed = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_MUNCOMPRESSED, false);
	//faidx
	//final SettingsModelString infasta = new SettingsModelString(SamToolsNodeModel.CFGKEY_INFASTA, "");
	//phase
	final SettingsModelIntegerBounded blocklength = new SettingsModelIntegerBounded(SamToolsNodeModel.CFGKEY_BLOCKLENGTH, 13, 1, Integer.MAX_VALUE);
	final SettingsModelString prefix = new SettingsModelString(SamToolsNodeModel.CFGKEY_PREFIX,"phase");
	final SettingsModelIntegerBounded hetphred = new SettingsModelIntegerBounded(SamToolsNodeModel.CFGKEY_HETPHRED, 37, 1, Integer.MAX_VALUE);
	final SettingsModelIntegerBounded minqual = new SettingsModelIntegerBounded(SamToolsNodeModel.CFGKEY_MINQUAL, 13, 1, Integer.MAX_VALUE);
	final SettingsModelIntegerBounded maxdepth = new SettingsModelIntegerBounded(SamToolsNodeModel.CFGKEY_MAXDEPTH, 256, 1, Integer.MAX_VALUE);
	final SettingsModelBoolean fixchimeras = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_FIXCHIMERAS, false);
	//final SettingsModelBoolean dropambig = new SettingsModelBoolean(SamToolsNodeModel.CFGKEY_DROPAMBIG, false);
	
	
	
    protected SamToolsNodeDialog() {
        
    	addPrefPageSetting(samtools, IBISKNIMENodesPlugin.SAMTOOLS);
    	addPrefPageSetting(refseqfile, IBISKNIMENodesPlugin.REF_GENOME);
    	
    	createNewGroup("Select utility");
    	addDialogComponent(new DialogComponentStringSelection(utility,"Select Utility", "cat", "faidx", "calmd", "fixmate", "flagstat","idxstats","merge", "phase", "reheader", "rmdup"));

    	
    	createNewGroup("SamTools");
    	addDialogComponent(new DialogComponentFileChooser(samtools, "SamTools", 0,""));
    	createNewGroup("Input (sorted) BAM file");
    	addDialogComponent(new DialogComponentFileChooser(bamfile, "Input BAM file", 0,".bam", ".BAM"));
    	createNewGroup("Reference sequence fasta file");
    	addDialogComponent(new DialogComponentFileChooser(refseqfile, "Reference sequence file", 0,""));
    	
    	createNewTab("cat");
       	setHorizontalPlacement(true);
    	createNewGroup("Header SAM file:");
    	addDialogComponent(new DialogComponentBoolean(useHeaderSAM, "Use header SAM"));
    	addDialogComponent(new DialogComponentFileChooser(headerSAM, "his_id0_samtools", 0,".sam", ".SAM"));
    	setHorizontalPlacement(false);
    	createNewGroup("Second input BAM file:");
    	addDialogComponent(new DialogComponentFileChooser(inBAM1, "Second BAM file", 0,".bam", ".BAM"));
      	//addDialogComponent(new DialogComponentFileChooser(inBAM2, "Second BAM file", 0,".bam", ".BAM"));
      	
      	//createNewTab("faidx");
      	//createNewGroup("Fasta file to index:");
      	//addDialogComponent(new DialogComponentFileChooser(infasta, "", 0, false));
      	
    	createNewTab("calmd");
    	addDialogComponent(new DialogComponentBoolean(changeIdentBases, "Change ident. bases to '='"));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(useCompression, "Output as BAM file"));
    	addDialogComponent(new DialogComponentStringSelection(compression, "","compressed", "uncompressed"));//use uncompressed for pipelining!
       	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentBoolean(inputIsSAM, "Input is SAM with header"));
    	addDialogComponent(new DialogComponentBoolean(modifyQual, "Modify the quality string"));
    	addDialogComponent(new DialogComponentBoolean(bqTag, "Compute BQ tag if qual. String modif., else cap baseQ by BAQ"));
       	addDialogComponent(new DialogComponentBoolean(extendedBAQ, "Extended BAQ for better sensitivity but lower specificity"));
    	/*setHorizontalPlacement(true);
       	addDialogComponent(new DialogComponentBoolean(doCapMapQual, "Use coefficient to cap mapping quality of poorly mapped reads"));
    	addDialogComponent(new DialogComponentNumber(capMapQual,"", 1));
       	setHorizontalPlacement(false);*/
       	
       	createNewTab("merge");
    	createNewGroup("Second input BAM file:");
    	addDialogComponent(new DialogComponentFileChooser(minbam1, "his_id1_samtools", 0,".bam", ".BAM"));
      	//addDialogComponent(new DialogComponentFileChooser(minbam2, "Second BAM file", 0,".bam", ".BAM"));
      	createNewGroup("Parameters:");
    	setHorizontalPlacement(true);
     	addDialogComponent(new DialogComponentBoolean(usemhfile, "Use lines of SAM file as header"));
     	addDialogComponent(new DialogComponentFileChooser(mhfile, "his_id2_samtools", 0,".sam", ".SAM"));
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
    	addDialogComponent(new DialogComponentNumber(blocklength,"Block length", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(hetphred,"Min. het phred-LOD", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(minqual,"Min. base quality in het calling", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(maxdepth,"Max. read depth", /*step*/ 1));
      	addDialogComponent(new DialogComponentBoolean(fixchimeras, "Do not attempt to fix chimeric reads"));
      	//addDialogComponent(new DialogComponentBoolean(dropambig, "Drop reads with ambiguous phase"));
 
      	
      	createNewTab("reheader");
      	createNewGroup("Header SAM file:");
    	addDialogComponent(new DialogComponentFileChooser(rehInSAM, "Header SAM file", 0,".sam", ".SAM"));
    	//createNewGroup("BAM file:");
    	//addDialogComponent(new DialogComponentFileChooser(rehInBAM, "BAM file", 0,".bam", ".BAM"));
    	
       	createNewTab("rmdup");
       	addDialogComponent(new DialogComponentBoolean(removeDup, "Remove duplicate for single-end reads"));
       	addDialogComponent(new DialogComponentBoolean(treatPE, "Treat paired-end reads and single-end reads."));
       	

      	
       	/**Checkboxes for the optional arguments**/
       
       	
       	
    	//checkboxes for calmd tab
    	utility.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				setEnabled(true, "cat");
				//calmd
				if(utility.getStringValue().equals("calmd")){
					changeIdentBases.setEnabled(true);
					compression.setEnabled(true);
					useCompression.setEnabled(true);
					inputIsSAM.setEnabled(true);
					modifyQual.setEnabled(true);
					bqTag.setEnabled(true);
					extendedBAQ.setEnabled(true);
					//doCapMapQual.setEnabled(true);
					setEnabled(true, "calmd");

					if(SamToolsNodeModel.useFastafile){
						refseqfile.setEnabled(true);
					}
			    	if(SamToolsNodeModel.useBamfile){
			    		bamfile.setEnabled(true);
			    	}
					
				}
				else{
					useCompression.setBooleanValue(false);
					//doCapMapQual.setBooleanValue(false);
					changeIdentBases.setEnabled(false);
					compression.setEnabled(false);
					useCompression.setEnabled(false);
					inputIsSAM.setEnabled(false);
					inputIsSAM.setBooleanValue(false);
					modifyQual.setEnabled(false);
					bqTag.setEnabled(false);
					extendedBAQ.setEnabled(false);
					//capMapQual.setEnabled(false);
					//doCapMapQual.setEnabled(false);
					setEnabled(false, "calmd");

				}
				//rmdup
				if(utility.getStringValue().equals("rmdup")){
					removeDup.setEnabled(true);
					treatPE.setEnabled(true);
					setEnabled(true, "rmdup");
					
			    	if(!SamToolsNodeModel.useSamtools){
						samtools.setEnabled(false);
					}
			    	if(SamToolsNodeModel.useBamfile){
			    		bamfile.setEnabled(true);
			    	}

				}
				else{
					removeDup.setEnabled(false);
					treatPE.setEnabled(false);
					setEnabled(false, "rmdup");
				}
				//cat
				if(utility.getStringValue().equals("cat")){
					useHeaderSAM.setEnabled(true);
					inBAM1.setEnabled(true);
					//inBAM2.setEnabled(true);
					setEnabled(true, "cat");
					
					refseqfile.setEnabled(false);
			    	if(SamToolsNodeModel.useBamfile){
			    		bamfile.setEnabled(true);
			    	}

				}
				else{
			    	inBAM1.setEnabled(false);
			    	useHeaderSAM.setEnabled(false);
			    	useHeaderSAM.setBooleanValue(false);
			    	headerSAM.setEnabled(false);
					setEnabled(false,"cat");
				}
				//faidx
				if(utility.getStringValue().equals("faidx")){
					if(SamToolsNodeModel.useFastafile){
						refseqfile.setEnabled(true);
					}
					if(SamToolsNodeModel.useBamfile){
						bamfile.setEnabled(false);
					}
					
				}
				
				//fixmate
				if(utility.getStringValue().equals("fixmate") || utility.getStringValue().equals("flagstat") || utility.getStringValue().equals("idxstats") || utility.getStringValue().equals("rmdup")){

					refseqfile.setEnabled(false);
			    	if(SamToolsNodeModel.useBamfile){
			    		bamfile.setEnabled(true);
			    	}

				}

				
				else{
					//useHeaderSAM.setBooleanValue(false);
					//useHeaderSAM.setEnabled(false);
					//headerSAM.setEnabled(false);
					//inBAM1.setEnabled(false);
					//inBAM2.setEnabled(false);
					setEnabled(false, "rmdup");
				}
				//reheader
				if(utility.getStringValue().equals("reheader")){
					rehInSAM.setEnabled(true);
					//rehInBAM.setEnabled(true);
					setEnabled(true, "reheader");
					
					refseqfile.setEnabled(false);
			    	if(SamToolsNodeModel.useBamfile){
			    		bamfile.setEnabled(true);
			    	}
					
				}
				else{
					rehInSAM.setEnabled(false);
					//rehInBAM.setEnabled(false);
					setEnabled(false, "reheader");
				}
				//merge
				if(utility.getStringValue().equals("merge")){
					mcompression.setEnabled(true);
					mforce.setEnabled(true);
					usemhfile.setEnabled(true);
					msorted.setEnabled(true);
					usemregion.setEnabled(true);
					mrgtag.setEnabled(true);
					muncompressed.setEnabled(true);
					minbam1.setEnabled(true);
					//minbam2.setEnabled(true);
					setEnabled(true, "merge");
					
					refseqfile.setEnabled(false);
			    	if(SamToolsNodeModel.useBamfile){
			    		bamfile.setEnabled(true);
			    	}

				}
				else{
					usemhfile.setBooleanValue(false);
					usemregion.setBooleanValue(false);
					mcompression.setEnabled(false);
					mforce.setEnabled(false);
					usemhfile.setEnabled(false);
					msorted.setEnabled(false);
					usemregion.setEnabled(false);
					mrgtag.setEnabled(false);
					muncompressed.setEnabled(false);
					minbam1.setEnabled(false);
					//minbam2.setEnabled(false);
					setEnabled(false, "merge");
				}

				//phase
				if(utility.getStringValue().equals("phase")){
			    	blocklength.setEnabled(true);
			    	prefix.setEnabled(true);
			    	hetphred.setEnabled(true);
			    	minqual.setEnabled(true);
			    	maxdepth.setEnabled(true);
			    	fixchimeras.setEnabled(true);
			    	//dropambig.setEnabled(true);
					setEnabled(true, "phase");
					
					refseqfile.setEnabled(false);
			    	if(SamToolsNodeModel.useBamfile){
			    		bamfile.setEnabled(true);
			    	}
					
				}
				else{
			    	blocklength.setEnabled(false);
			    	prefix.setEnabled(false);
			    	hetphred.setEnabled(false);
			    	minqual.setEnabled(false);
			    	maxdepth.setEnabled(false);
			    	fixchimeras.setEnabled(false);
			    	//dropambig.setEnabled(false);
					setEnabled(false, "phase");
				}
				
			}
		});
    	/*doCapMapQual.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(doCapMapQual.getBooleanValue()){
					capMapQual.setEnabled(true);
				}
				else{
					capMapQual.setEnabled(false);
				}
			}
		});*/
    	useCompression.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(useCompression.getBooleanValue()){
					compression.setEnabled(true);
				}
				else{
					compression.setEnabled(false);
				}
			}
		});
    	useHeaderSAM.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(useHeaderSAM.getBooleanValue()){
					headerSAM.setEnabled(true);
				}
				else{
					headerSAM.setEnabled(false);
				}
			}
		});
    	usemhfile.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usemhfile.getBooleanValue()){
					mhfile.setEnabled(true);
				}
				else{
					mhfile.setEnabled(false);
				}
			}
		});
    	usemregion.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usemregion.getBooleanValue()){
					mregion.setEnabled(true);
				}
				else{
					mregion.setEnabled(false);
				}
			}
		});
    	
    	



    	
    }

//	@Override
//	protected void updatePrefs() {
//		if(usePrefPage.getBooleanValue()) {
//			String toolPath = IBISKNIMENodesPlugin.getDefault().getToolPathPreference("samtools");
//		    if(toolPath != null && !toolPath.equals("")) {
//		    	samtools.setStringValue(toolPath);
//		    	samtools.setEnabled(false);
//		    } else {
//		    	samtools.setEnabled(true);
//		    }
//		    
//	    	String refGenome = IBISKNIMENodesPlugin.getDefault().getRefGenomePreference();
//	    	if(refGenome != null && !refGenome.equals("")) {
//	    		refseqfile.setStringValue(refGenome);
//	    		refseqfile.setEnabled(false);
//	    	} else {
//	    		refseqfile.setEnabled(true);
//	    	}
//		    
//		} else {
//			samtools.setEnabled(true);
//			refseqfile.setEnabled(true);
//		}
//	}

}

