package de.helmholtz_muenchen.ibis.ngs.bfast;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDouble;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;



/**
 * <code>NodeDialog</code> for the "RunBfast" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Sebastian
 */
public class BfastNodeDialog extends DefaultNodeSettingsPane {

	final SettingsModelBoolean usekeysize = new SettingsModelBoolean(BfastNodeModel.CFGKEY_USEKEYSIZE, false);
	final SettingsModelIntegerBounded keysize = new SettingsModelIntegerBounded(BfastNodeModel.CFGKEY_KEYSIZE, 4, 1, Integer.MAX_VALUE);
	//indexing
	final SettingsModelBoolean usesplitdepth = new SettingsModelBoolean(BfastNodeModel.CFGKEY_USESPLITDEPTH, false);
	final SettingsModelIntegerBounded splitdepth = new SettingsModelIntegerBounded(BfastNodeModel.CFGKEY_SPLITDEPTH, 0, 0, Integer.MAX_VALUE);
	final SettingsModelBoolean usetmpdir = new SettingsModelBoolean(BfastNodeModel.CFGKEY_USETMPDIR, false);
	final SettingsModelString tmpdir = new SettingsModelString(
			BfastNodeModel.CFGKEY_TMPDIR,"");
	
	final SettingsModelIntegerBounded startcontig = new SettingsModelIntegerBounded(
			BfastNodeModel.CFGKEY_STARTCONTIG,1,1,Integer.MAX_VALUE);
	final SettingsModelIntegerBounded endcontig = new SettingsModelIntegerBounded(
			BfastNodeModel.CFGKEY_ENDCONTIG,1,1,Integer.MAX_VALUE);
	final SettingsModelIntegerBounded startpos = new SettingsModelIntegerBounded(
			BfastNodeModel.CFGKEY_STARTPOS,1,1,Integer.MAX_VALUE);
	final SettingsModelIntegerBounded endpos = new SettingsModelIntegerBounded(
			BfastNodeModel.CFGKEY_ENDPOS,1,1,Integer.MAX_VALUE);
	final SettingsModelBoolean repeatmasker = new SettingsModelBoolean(BfastNodeModel.CFGKEY_REPEATMASKER, false);
	final SettingsModelString exonsfile = new SettingsModelString(BfastNodeModel.CFGKEY_EXONSFILE,"");
	final SettingsModelBoolean usecontigs = new SettingsModelBoolean(BfastNodeModel.CFGKEY_USECONTIGS, false);
	final SettingsModelBoolean usepos = new SettingsModelBoolean(BfastNodeModel.CFGKEY_USEPOS, false);
	final SettingsModelBoolean useexonsfile = new SettingsModelBoolean(BfastNodeModel.CFGKEY_USEEXONSFILE, false);
	
	//CALs
	final SettingsModelBoolean usemainindexes = new SettingsModelBoolean(BfastNodeModel.CFGKEY_USEMAININDEXES, false);
	final SettingsModelBoolean usesecindexes = new SettingsModelBoolean(BfastNodeModel.CFGKEY_USESECINDEXES, false);
	final SettingsModelBoolean useoffsets = new SettingsModelBoolean(BfastNodeModel.CFGKEY_USEOFFSETS, false);
	final SettingsModelString mainindexes = new SettingsModelString(
			BfastNodeModel.CFGKEY_MAININDEXES,"");
	final SettingsModelString secindexes = new SettingsModelString(
			BfastNodeModel.CFGKEY_SECINDEXES,"");
	final SettingsModelString offsets = new SettingsModelString(
			BfastNodeModel.CFGKEY_OFFSETS,"");
	final SettingsModelBoolean usemaxreads = new SettingsModelBoolean(BfastNodeModel.CFGKEY_USEMAXREADS, false);
	final SettingsModelIntegerBounded maxreads = new SettingsModelIntegerBounded(BfastNodeModel.CFGKEY_MAXREADS, 500000, 1, Integer.MAX_VALUE);	
	
	//Alignment&Postprocessing
	final SettingsModelString gap = new SettingsModelString(BfastNodeModel.CFGKEY_GAPPEDALIGN,"gapped");
	final SettingsModelIntegerBounded avg = new SettingsModelIntegerBounded(BfastNodeModel.CFGKEY_AVGMISQUAL, BfastNodeModel.DEFAULT_AVGMISQUAL, 0, Integer.MAX_VALUE);
	
	final SettingsModelBoolean useavgdev = new SettingsModelBoolean(BfastNodeModel.CFGKEY_USEAVGDEV, false);
	final SettingsModelIntegerBounded inssizeavg = new SettingsModelIntegerBounded(
				BfastNodeModel.CFGKEY_INSSIZEAVG,
				BfastNodeModel.DEFAULT_INSSIZEAVG,
				0,Integer.MAX_VALUE);
	final SettingsModelIntegerBounded insstddev = new SettingsModelIntegerBounded(
				BfastNodeModel.CFGKEY_INSSTDDEV,
				BfastNodeModel.DEFAULT_INSSTDDEV,
				1,Integer.MAX_VALUE);

	 
    /**
     * New pane for configuring the RunBfast node.
     */
    protected BfastNodeDialog() {
    	//setDefaultTabTitle("Ref.Genome");//Occurs at wrong position?
    	createNewGroup("Bfast");
    	addDialogComponent(new DialogComponentFileChooser(new SettingsModelString(
    			BfastNodeModel.CFGKEY_BFASTEXE,null), "Enter path to bfast exec.", 0, false));
    	
    	createNewGroup("Input fasta file (reference sequence)");
    	addDialogComponent(new DialogComponentFileChooser(new SettingsModelString(
    			BfastNodeModel.CFGKEY_INPUTFASTA,null), "Enter path to input .fasta", 0, false));
    	
    	createNewGroup("Ref. genome model type");
    	addDialogComponent(new DialogComponentStringSelection(
    			new SettingsModelString(BfastNodeModel.CFGKEY_MODELTYPE,"Nucleotide space"), "Model Type", "Nucleotide space", "Color space"));
    	
    	createNewGroup("Threading");
    	addDialogComponent(new DialogComponentNumber(
        		new SettingsModelIntegerBounded(BfastNodeModel.CFGKEY_THREADS, 1,1,Integer.MAX_VALUE),
        					"Number of Threads to use:", /*step*/ 1));

    	createNewTab("Indexing");
    	createNewGroup("Mask");
    	addDialogComponent(new DialogComponentString(new SettingsModelString(BfastNodeModel.CFGKEY_MASK,""), "Mask"));
    	createNewGroup("Hashwidth");
    	addDialogComponent(new DialogComponentNumber(
    		new SettingsModelIntegerBounded(BfastNodeModel.CFGKEY_HASHWIDTH, 12,1,Integer.MAX_VALUE),
    					"Hashwidth:", /*step*/ 1));
    	
    	createNewGroup("Split depth");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usesplitdepth, "Give split index depth (4^x)"));
    	addDialogComponent(new DialogComponentNumber(splitdepth, "Split depth", /*step*/ 1));
    	setHorizontalPlacement(false);
       	createNewGroup("Temporary directory");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usetmpdir, "Specify a temp. directory"));
      	addDialogComponent(new DialogComponentFileChooser(tmpdir, "Choose temp. directory.", 0, true));
    	setHorizontalPlacement(false);
    	
       	createNewGroup("Mask repeats");
    	addDialogComponent(new DialogComponentBoolean(repeatmasker, "Ignore lower case bases when creating the indexes"));
    	
    	createNewGroup("Start/End Contigs for indexing and first pos. in first contig, last pos. in last contig");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usecontigs, "Give start/end contigs for indexing"));
    	addDialogComponent(new DialogComponentNumber(startcontig, "startcontig", /*step*/ 1));
       	addDialogComponent(new DialogComponentNumber(endcontig, "endcontig", /*step*/ 1));
    	setHorizontalPlacement(false);
    	
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usepos, "Specify first/last position"));
    	addDialogComponent(new DialogComponentNumber(startpos, "startpos.", /*step*/ 1));
       	addDialogComponent(new DialogComponentNumber(endpos, "endpos.", /*step*/ 1));
    	setHorizontalPlacement(false);
    	
       	createNewGroup("Specify exon-like ranges to include in index");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(useexonsfile, "Use exons file"));
      	addDialogComponent(new DialogComponentFileChooser(exonsfile, "", 0, false));
    	setHorizontalPlacement(false);


    	
    	//New Tab
    	createNewTab("CALs");
    	//createNewGroup("Input file containing reads:");
    	//addDialogComponent(new DialogComponentFileChooser(new SettingsModelString(	RunBfastNodeModel.CFGKEY_READSINPUTFILE,null), "Enter path to file containing reads", 0, false));
    	//createNewGroup("Further options:");
    	createNewGroup("Memory policy");
    	addDialogComponent(new DialogComponentStringSelection(
    			new SettingsModelString(BfastNodeModel.CFGKEY_LOADINDEXES,"no"), "Load indexes in memory","yes", "no"));
    	createNewGroup("Compression");
    	addDialogComponent(new DialogComponentStringSelection(
    			new SettingsModelString(BfastNodeModel.CFGKEY_COMPRESSTYPE,"none"), "Reads compression type","none", "bz2", "gz"));
    	createNewGroup("Strand selection");
    	addDialogComponent(new DialogComponentStringSelection(
    			new SettingsModelString(BfastNodeModel.CFGKEY_STRAND,"Both"), "Search which strand","both", "forward", "reverse"));
    	
    	createNewGroup("Main indexes");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usemainindexes, "Main Indexes (comma sep.)"));
    	addDialogComponent(new DialogComponentString(mainindexes, ""));
       	setHorizontalPlacement(false);
    	createNewGroup("Secondary indexes");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usesecindexes, "Secondary Indexes (comma sep.)"));
    	addDialogComponent(new DialogComponentString(secindexes, ""));
       	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	createNewGroup("Offsets");
    	addDialogComponent(new DialogComponentBoolean(useoffsets, "Offsets (range or comma sep.)"));
    	addDialogComponent(new DialogComponentString(offsets, ""));
       	setHorizontalPlacement(false);
    	
    	createNewGroup("Read restriction");
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    				BfastNodeModel.CFGKEY_STARTREADNUM, 1, 1, Integer.MAX_VALUE),
    				"Start Read Number:", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    				BfastNodeModel.CFGKEY_ENDREADNUM, 2147483647, 1, Integer.MAX_VALUE),
    				"End Read Number:", /*step*/ 1));
    	createNewGroup("Keysize");
       	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usekeysize, "Keysize:"));
    	addDialogComponent(new DialogComponentNumber(keysize,"", /*step*/ 1));
       	setHorizontalPlacement(false);
    	
    	createNewGroup("Maximal keymatches");
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    				BfastNodeModel.CFGKEY_MAXKEYMATCHES, 8, 1, Integer.MAX_VALUE),
    				"Max. Keymatches:", /*step*/ 1));
    	createNewGroup("Keymiss fraction");
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelDouble(BfastNodeModel.CFGKEY_KEYMISSFRACTION, 1.0),
    				"Keymiss fraction:", /*step*/ 0.01));
    	createNewGroup("Max. number of CALs per read");
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    				BfastNodeModel.CFGKEY_MAXNUMMATCHES, 384, 1, Integer.MAX_VALUE),
    				"Max. matches:", /*step*/ 1));
    	createNewGroup("Number of reads to load at a time");
       	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usemaxreads, "Limit number of reads to load at a time"));
    	addDialogComponent(new DialogComponentNumber(maxreads,"", /*step*/ 250000));
       	setHorizontalPlacement(false);

		/**************************************
		*		Alignment&postprocessing	  *
		**************************************/
    	createNewTab("Align");
    	createNewGroup("Type");
       	addDialogComponent(new DialogComponentStringSelection(gap,"Local Alignment", "gapped","ungapped"));
    	createNewGroup("Mask constraints");
       	addDialogComponent(new DialogComponentStringSelection(
    			new SettingsModelString(BfastNodeModel.CFGKEY_MASKCONSTRAINTS,""),
    			"Use mask constraints from the match step", "Yes","No"));
    	/**createNewGroup("Read start/stop");
       	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					BfastNodeModel.CFGKEY_READSTART, 
    					BfastNodeModel.DEFAULT_READSTART, 0, Integer.MAX_VALUE),
    					"Specify the read to begin with:",  1));
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					BfastNodeModel.CFGKEY_READSTOP, 
    					BfastNodeModel.DEFAULT_READSTOP, 0, Integer.MAX_VALUE),
    					"Specify the read to stop with:",  1));**/
    	createNewGroup("Offset");
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					BfastNodeModel.CFGKEY_OFFSET, 
    					BfastNodeModel.DEFAULT_OFFSET, 0, Integer.MAX_VALUE),
    					"Specify the number of bases before and after\n the match to include in the reference genome", /*step*/ 1));
    	createNewGroup("Max. candidates");
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					BfastNodeModel.CFGKEY_MAXMATCHES, 
    					BfastNodeModel.DEFAULT_MAXMATCHES, 0, Integer.MAX_VALUE),
    					"Specify the maximum number of candidates to initiate alignment for a given match: ", /*step*/ 1));
    	createNewGroup("Mismatch quality");
    	addDialogComponent(new DialogComponentNumber(
    			avg,
    					"Specify the average mismatch quality:", /*step*/ 1));
    	
    	createNewTab("Postprocessing");
    	createNewGroup("Algorithm");
    	addDialogComponent(new DialogComponentStringSelection(
    			new SettingsModelString(BfastNodeModel.CFGKEY_ALGO,BfastNodeModel.DEFAULT_ALGO), "" +
    					" Specify the algorithm to choose the alignment for each end of the read",
    					"No filtering", "All alignments","Uniquely aligned reads only", "Unique alignment with best score","All alignments with best score"));
    	createNewGroup("Pairing");
    	addDialogComponent(new DialogComponentStringSelection(
    			new SettingsModelString(BfastNodeModel.CFGKEY_PAIRING,BfastNodeModel.DEFAULT_PAIRING), "" +
    					"Specify the pairing orientation:",
    					"paired ends", "mate pairs", "no pairing"));
    	createNewGroup("Minimal map quality");
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					BfastNodeModel.CFGKEY_MINMAPQUAL, 
    					BfastNodeModel.DEFAULT_MINMAPQUAL, Integer.MIN_VALUE, Integer.MAX_VALUE),
    					"Specify to remove low mapping quality alignments:", /*step*/ 1));
    	createNewGroup("Minimal score");
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    					BfastNodeModel.CFGKEY_MINNORMSCORE, 
    					BfastNodeModel.DEFAULT_MINNORMSCORE, Integer.MIN_VALUE, Integer.MAX_VALUE),
    					"Specify to remove low (alignment) scoring alignments:", /*step*/ 1));
    	createNewGroup("Mean insert size & standard deviation");
       	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(useavgdev, "Specify insert size:"));
    	addDialogComponent(new DialogComponentNumber(inssizeavg,
    					"Specify the mean insert size to use when pairing:", /*step*/ 1));
       	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentNumber(insstddev,"Specify the standard deviation of the insert size to use when pairing:", /*step*/ 1));

    	/** Event Listeners **/
    	useavgdev.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(useavgdev.getBooleanValue()){
					insstddev.setEnabled(true);
					inssizeavg.setEnabled(true);
				}
				else{
					insstddev.setEnabled(false);
					inssizeavg.setEnabled(false);
				}
			}
		});
    	usemainindexes.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usemainindexes.getBooleanValue()){
					mainindexes.setEnabled(true);
				}
				else{
					mainindexes.setEnabled(false);
				}
			}
		});
    	usesecindexes.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usesecindexes.getBooleanValue()){
					secindexes.setEnabled(true);
				}
				else{
					secindexes.setEnabled(false);
				}
			}
		});
    	useoffsets.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(useoffsets.getBooleanValue()){
					offsets.setEnabled(true);
				}
				else{
					offsets.setEnabled(false);
				}
			}
		});
    	usekeysize.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usekeysize.getBooleanValue()){
					keysize.setEnabled(true);
				}
				else{
					keysize.setEnabled(false);
				}
			}
		});
    	usesplitdepth.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usesplitdepth.getBooleanValue()){
					splitdepth.setEnabled(true);
				}
				else{
					splitdepth.setEnabled(false);
				}
			}
		});
    	usetmpdir.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usetmpdir.getBooleanValue()){
					tmpdir.setEnabled(true);
				}
				else{
					tmpdir.setEnabled(false);
				}
			}
		});
    	usecontigs.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usecontigs.getBooleanValue()){
					startcontig.setEnabled(true);
					endcontig.setEnabled(true);
					useexonsfile.setBooleanValue(false);
					exonsfile.setEnabled(false);
				}
				else{
					startcontig.setEnabled(false);
					endcontig.setEnabled(false);
				}
			}
		});
    	usepos.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usepos.getBooleanValue()){
					startpos.setEnabled(true);
					endpos.setEnabled(true);
					useexonsfile.setBooleanValue(false);
					exonsfile.setEnabled(false);
				}
				else{
					startpos.setEnabled(false);
					endpos.setEnabled(false);
				}
			}
		});
    	useexonsfile.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(useexonsfile.getBooleanValue()){
					usepos.setBooleanValue(false);
					usecontigs.setBooleanValue(false);
					startpos.setEnabled(false);
					endpos.setEnabled(false);
					startcontig.setEnabled(false);
					endcontig.setEnabled(false);
					exonsfile.setEnabled(true);
				}
				else{
					exonsfile.setEnabled(false);
				}
			}
		});
    	usemaxreads.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usemaxreads.getBooleanValue()){
					maxreads.setEnabled(true);

				}
				else{
					maxreads.setEnabled(false);
				}
			}
		});
    	

    	
    }
}