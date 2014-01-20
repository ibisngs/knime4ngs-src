package de.helmholtz_muenchen.ibis.ngs.annovar;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "Annovar" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class AnnovarNodeDialog extends DefaultNodeSettingsPane {

	/**General**/
	final SettingsModelString queryfile = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_QUERYFILE,"");
	final SettingsModelString databaselocation = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_DATABASELOCATION,"");
	final SettingsModelString path2annovar = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_PATH2ANNOVAR,"");
	final SettingsModelString method = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_METHOD,"");
	final SettingsModelBoolean usetablename = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_USETABLENAME, false);
	final SettingsModelString tablename = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_TABLENAME,"");
	
	/**Arguments to control in-/output**/
    /*--outfile <file>            output file prefix
    --zerostart                 input query file uses half-open zero-start coordinate
    --dbtype <string>           database type
    --buildver <string>         genome build version (default: hg18 for human)
    --gff3dbfile <file>         specify the GFF3 DB file used in region-based annotation
    --genericdbfile <file>      specify the generic DB file used in filter-based annotation
    --vcfdbfile <file>          specify the DB file in VCF format in filter-based annotation
    --bedfile <file>            specify a BED file in region-based annotation
    --time                      print out local time during program run
    --separate                  separately print out all function of a variant (default: one line per variant)
    --colsWanted <string>       specify which columns to output in -regionanno by comma-delimited numbers
    --comment                   print out comment line (those starting with #) in output files 
    --scorecolumn <int>         the column with scores in database file (for region-based annotation)
    --exonsort                  sort the exon number in output line (for gene-based annotation)
    --transcript_function       use transcript name rather than gene name in gene-based annotation output
    --hgvs                      use HGVS format for exonic annotation (c.122C>T rather than c.C122T)
    --otherinfo                 in filter-based annotation, print out additional columns in database file
    --infoasscore             in filter-based annotation, use INFO field in VCF file as score in output
    --seq_padding               if set, create a new file with cDNA sequence padded by this much either side*/
	//final SettingsModelBoolean useoutfile = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_USEOUTFILE, false);
	final SettingsModelString outfile = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_OUTFILE,AnnovarNodeModel.DEFAULT_PRECEDENCE);
	//final SettingsModelBoolean usedbtype = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_USEDBTYPE, false);
	final SettingsModelString dbtype = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_DBTYPE,"");
	//final SettingsModelBoolean usebuildver = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_USEBUILDVER, false);
	final SettingsModelString buildver = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_BUILDVER,"");
	final SettingsModelBoolean usegff3dbfile = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_USEGFF3DBFILE, false);
	final SettingsModelString gff3dbfile = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_GFF3DBFILE,"");
	final SettingsModelBoolean usegenericdbfile = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_USEGENERICDBFILE, false);
	final SettingsModelString genericdbfile = new SettingsModelString(AnnovarNodeModel.CFGKEY_GENERICDBFILE,"");
	final SettingsModelBoolean usevcfdbfile = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_USEVCFDBFILE, false);
	final SettingsModelString vcfdbfile = new SettingsModelString(AnnovarNodeModel.CFGKEY_VCFDBFILE,"");
	final SettingsModelBoolean usebedfile = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_USEBEDFILE, false);
	final SettingsModelString bedfile = new SettingsModelString(AnnovarNodeModel.CFGKEY_BEDFILE,"");
	final SettingsModelBoolean separate = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_SEPARATE, false);
	final SettingsModelBoolean usecolswanted = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_USECOLSWANTED, false);
	final SettingsModelString colswanted = new SettingsModelString(AnnovarNodeModel.CFGKEY_COLSWANTED,"");
	final SettingsModelBoolean comment = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_COMMENT, false);
	final SettingsModelBoolean usescorecolumn = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_USESCORECOLUMN, false);
	final SettingsModelIntegerBounded scorecolumn = new SettingsModelIntegerBounded(
			AnnovarNodeModel.CFGKEY_SCORECOLUMN,0,0,Integer.MAX_VALUE);
	final SettingsModelBoolean exonsort = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_EXONSORT, false);
	final SettingsModelBoolean transcript_function = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_TRANSCRIPTFUNCTION, false);
	final SettingsModelBoolean hgvs = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_HGVS, false);
	final SettingsModelBoolean seq_padding = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_SEQPADDING, false);
	//filterbased
	final SettingsModelBoolean otherinfo = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_OTHERINFO, false);
	final SettingsModelBoolean infoasscore = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_INFOASSCORE, false);

	/**Arguments to fine-tune the annotation procedure**/
	
	  /**  Arguments to fine-tune the annotation procedure
    --batchsize <int>           batch size for processing variants per batch (default: 5m)
    --genomebinsize <int>       bin size to speed up search (default: 100k for -geneanno, 10k for -regionanno)
    --expandbin <int>           check nearby bin to find neighboring genes (default: 2m/genomebinsize)
    !--neargene <int>            distance threshold to define upstream/downstream of a gene
    !--score_threshold <float>   minimum score of DB regions to use in annotation
    --reverse                   reverse directionality to compare to score_threshold
   ! --normscore_threshold <float> minimum normalized score of DB regions to use in annotation
    --rawscore                  output includes the raw score (not normalized score) in UCSC Browser Track
    --minqueryfrac <float>      minimum percentage of query overlap to define match to DB (default: 0)
    !--splicing_threshold <int>  distance between splicing variants and exon/intron boundary (default: 2)
    !--indel_splicing_threshold <int>    if set, use this value for allowed indel size for splicing variants (default: --splicing_threshold)
    --maf_threshold <float>     filter 1000G variants with MAF above this threshold (default: 0)
    --sift_threshold <float>    SIFT threshold for deleterious prediction (default: 0.05)
    --precedence <string>       comma-delimited to specify precedence of variant function (default: exonic>intronic...)
    --indexfilter_threshold <float>     controls whether filter-based annotation use index if this fraction of bins need to be scanned (default: 0.9)
**/
	final SettingsModelIntegerBounded batchsize = new SettingsModelIntegerBounded(
			AnnovarNodeModel.CFGKEY_BATCHSIZE,AnnovarNodeModel.DEFAULT_BATCHSIZE,1,Integer.MAX_VALUE);
	final SettingsModelIntegerBounded genomebinsize = new SettingsModelIntegerBounded(
			AnnovarNodeModel.CFGKEY_GENOMEBINSIZE,AnnovarNodeModel.DEFAULT_GENOMEBINSIZE,1,Integer.MAX_VALUE);
	final SettingsModelIntegerBounded expandbin = new SettingsModelIntegerBounded(
			AnnovarNodeModel.CFGKEY_EXPANDBIN,AnnovarNodeModel.DEFAULT_EXPANDBIN,1,Integer.MAX_VALUE);
	final SettingsModelBoolean useneargene = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_USENEARGENE, false);
	final SettingsModelIntegerBounded neargene = new SettingsModelIntegerBounded(
			AnnovarNodeModel.CFGKEY_NEARGENE,0,0,Integer.MAX_VALUE);
	final SettingsModelBoolean usescorethreshold = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_USESCORETHRESHOLD, false);
	final SettingsModelDoubleBounded score_threshold = new SettingsModelDoubleBounded(
			AnnovarNodeModel.CFGKEY_SCORETHRESHOLD,0,0,Double.MAX_VALUE);
	final SettingsModelBoolean reverse = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_REVERSE, false);
	final SettingsModelBoolean usenormscorethreshold = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_USENORMSCORETHRESHOLD, false);
	final SettingsModelIntegerBounded normscore_threshold = new SettingsModelIntegerBounded(
			AnnovarNodeModel.CFGKEY_NORMSCORETHRESHOLD,0,0,Integer.MAX_VALUE);
	final SettingsModelBoolean rawscore = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_RAWSCORE, false);
	final SettingsModelDoubleBounded minqueryfrac = new SettingsModelDoubleBounded(
			AnnovarNodeModel.CFGKEY_MINQUERYFRAC,AnnovarNodeModel.DEFAULT_MINQUERYFRAC,0,Double.MAX_VALUE);
	final SettingsModelBoolean usesplicingthreshold = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_USESPLICINGTHRESHOLD, false);
	final SettingsModelIntegerBounded splicing_threshold = new SettingsModelIntegerBounded(
			AnnovarNodeModel.CFGKEY_SPLICINGTHRESHOLD,AnnovarNodeModel.DEFAULT_SPLICINGTHRESHOLD,0,Integer.MAX_VALUE);
	final SettingsModelBoolean useindelsplicingthreshold = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_USEINDELSPLICINGTHRESHOLD, false);
	final SettingsModelIntegerBounded indel_splicing_threshold = new SettingsModelIntegerBounded(
			AnnovarNodeModel.CFGKEY_INDELSPLICINGTHRESHOLD,1,1,Integer.MAX_VALUE);
	final SettingsModelDoubleBounded maf_threshold = new SettingsModelDoubleBounded(
			AnnovarNodeModel.CFGKEY_MAFTHRESHOLD,AnnovarNodeModel.DEFAULT_MAFTHRESHOLD,0,Double.MAX_VALUE);
	final SettingsModelBoolean usesift_threshold = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_USESIFTTHRESHOLD, false);
	final SettingsModelDoubleBounded sift_threshold = new SettingsModelDoubleBounded(
			AnnovarNodeModel.CFGKEY_SIFTTHRESHOLD,AnnovarNodeModel.DEFAULT_SIFTTHRESHOLD,0,Double.MAX_VALUE);
	final SettingsModelString precedence = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_PRECEDENCE,AnnovarNodeModel.DEFAULT_PRECEDENCE);
	final SettingsModelDoubleBounded indexfilter_threshold = new SettingsModelDoubleBounded(
			AnnovarNodeModel.CFGKEY_INDEXFILTERTHRESHOLD,AnnovarNodeModel.DEFAULT_INDEXFILTERTHRESHOLD,0,Double.MAX_VALUE);
	
	/**arguments to control memory usage**/
	final SettingsModelIntegerBounded memfree = new SettingsModelIntegerBounded(
			AnnovarNodeModel.CFGKEY_MEMFREE,AnnovarNodeModel.DEFAULT_MEMFREE,100000,Integer.MAX_VALUE);
	final SettingsModelIntegerBounded memtotal = new SettingsModelIntegerBounded(
			AnnovarNodeModel.CFGKEY_MEMTOTAL,AnnovarNodeModel.DEFAULT_MEMTOTAL,0,Integer.MAX_VALUE);
	final SettingsModelBoolean usechromosome = new SettingsModelBoolean(AnnovarNodeModel.CFGKEY_USECHROMOSOME, false);
	final SettingsModelString chromosome = new SettingsModelString(
			AnnovarNodeModel.CFGKEY_CHROMOSOME,"");
	
    /**
     * New pane for configuring the Annovar node.
     */
    protected AnnovarNodeDialog() {
    	createNewGroup("Directory containing Annovar scripts");
    	addDialogComponent(new DialogComponentFileChooser(path2annovar, "his_ann_ID0", 0, true));
    	createNewGroup("Tool");
    	addDialogComponent(new DialogComponentStringSelection(method,"Select Utility", "geneanno", "regionanno" ,"filter"));
    	createNewGroup("Query file");
    	addDialogComponent(new DialogComponentFileChooser(queryfile, "his_ann_ID1", 0, false));
    	createNewGroup("Table name (instead of query file)");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usetablename, "Specify table name"));
    	addDialogComponent(new DialogComponentString(tablename, ""));
    	setHorizontalPlacement(false);
    	createNewGroup("Database location");
    	addDialogComponent(new DialogComponentFileChooser(databaselocation, "his_ann_ID2", 0, true));

    	createNewTab("In/Out");
    	createNewGroup("Specify output directory");
    	addDialogComponent(new DialogComponentFileChooser(outfile, "his_ann_ID3", 0, true));
    	createNewGroup("Database type");
    	//setHorizontalPlacement(true);
    	//addDialogComponent(new DialogComponentBoolean(usedbtype, "set database type"));
       	addDialogComponent(new DialogComponentString(dbtype, "specify database type"));
    	//setHorizontalPlacement(false);
    	createNewGroup("Build version");
    	//setHorizontalPlacement(true);
    	//addDialogComponent(new DialogComponentBoolean(usebuildver, "Genome build version"));
       	addDialogComponent(new DialogComponentString(buildver, "specify build version"));
    	//setHorizontalPlacement(false);
    	createNewGroup("Oher options");
    	addDialogComponent(new DialogComponentBoolean(separate, "Print variant functions separately"));
    	addDialogComponent(new DialogComponentBoolean(comment, "Print comment lines in output files"));

    	addDialogComponent(new DialogComponentBoolean(hgvs, "use HGVS format for exonic annotation"));
    	addDialogComponent(new DialogComponentBoolean(seq_padding, "create a new file with cDNA sequence padded by this much either side"));
    	createNewGroup("Gene-based options");
    	addDialogComponent(new DialogComponentBoolean(exonsort, "sort the exon number in output line"));
    	addDialogComponent(new DialogComponentBoolean(transcript_function, "use transcript name rather than gene name in output"));
    	
       //region based
    	createNewGroup("Region-based options");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usegff3dbfile, "Specify gff3 db file"));
    	addDialogComponent(new DialogComponentFileChooser(gff3dbfile, "", 0, false));
    	setHorizontalPlacement(false);
    	
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usebedfile, "Specify bedfile"));
    	addDialogComponent(new DialogComponentFileChooser(bedfile, "", 0, false));
    	setHorizontalPlacement(false);
    	
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usecolswanted, "specify which columns to output (comma-delimited num.)"));
       	addDialogComponent(new DialogComponentString(colswanted, ""));
    	setHorizontalPlacement(false);
    	
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usescorecolumn, "specify column with scores in db file"));
    	addDialogComponent(new DialogComponentNumber(scorecolumn,	"", /*step*/ 1));
    	setHorizontalPlacement(false);
    	
    	//filter based
    	createNewGroup("Filter-based options");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usegenericdbfile, "Specify generic db file"));
    	addDialogComponent(new DialogComponentFileChooser(genericdbfile, "", 0, false));
    	setHorizontalPlacement(false);
    	
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usevcfdbfile, "Specify vcf db file"));
    	addDialogComponent(new DialogComponentFileChooser(vcfdbfile, "", 0, false));
    	setHorizontalPlacement(false);

    	addDialogComponent(new DialogComponentBoolean(otherinfo, "print out additional columns in database file"));
    	addDialogComponent(new DialogComponentBoolean(infoasscore, "use INFO field in VCF file as score in output"));
    	/**Fine-tuning**/
    	createNewTab("Fine-tuning");
    	createNewGroup("Batch");
    	addDialogComponent(new DialogComponentNumber(batchsize, "Batch size", /*step*/ 1));
    	createNewGroup("Bins");
       	addDialogComponent(new DialogComponentNumber(genomebinsize, "Genome bin size", /*step*/ 1));
       	addDialogComponent(new DialogComponentNumber(expandbin, "Expand bin", /*step*/ 1));
    	createNewGroup("Distance threshold to define upstream/downstream of a gene");
       	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(useneargene, "use threshold"));
    	addDialogComponent(new DialogComponentNumber(neargene, "distance", /*step*/ 1));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentBoolean(reverse, "reverse directionality to compare to score threshold"));
     	createNewGroup("Rawscore output");
    	addDialogComponent(new DialogComponentBoolean(rawscore, "include raw score in output"));

    	createNewGroup("Minimum score of DB regions (not available for gene anno)");
       	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usescorethreshold, "use score threshold"));
    	addDialogComponent(new DialogComponentNumber(score_threshold, "set threshold", /*step*/ 0.01));
    	setHorizontalPlacement(false);
    	//gene anno
    	createNewGroup("Gene-based");
    	addDialogComponent(new DialogComponentString(precedence, "precedence"));
       	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usesplicingthreshold, "use splicing threshold"));
    	addDialogComponent(new DialogComponentNumber(splicing_threshold, "set threshold", /*step*/ 1));
    	setHorizontalPlacement(false);
       	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(useindelsplicingthreshold, "use indel splicing threshold"));
    	addDialogComponent(new DialogComponentNumber(indel_splicing_threshold, "set threshold", /*step*/ 1));
    	setHorizontalPlacement(false);
    	
    	//region anno
    	createNewGroup("Region-based");
       	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usenormscorethreshold, "minimum normalized score of DB regions"));
    	addDialogComponent(new DialogComponentNumber(normscore_threshold, "score threshold", /*step*/ 1));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentNumber(minqueryfrac, "Min. query frac.", /*step*/ 0.01));
    	//filter options
    	createNewGroup("Filter-based");
    	addDialogComponent(new DialogComponentNumber(maf_threshold, "MAF threshold", /*step*/ 0.01));
       	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usesift_threshold, "Set sift threshold (only if dbtype=avsift)"));
    	addDialogComponent(new DialogComponentNumber(sift_threshold, "SIFT threshold", /*step*/ 0.01));
       	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentNumber(indexfilter_threshold, "Indexfilter threshold", /*step*/ 0.01));
    	
    	createNewTab("Memory");
       	addDialogComponent(new DialogComponentNumber(memfree, "Minimal required memory [kb]", /*step*/ 1));
       	addDialogComponent(new DialogComponentNumber(memtotal, "Limit memory used by Annovar [kb]", /*step*/ 1));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usechromosome, "Examine specific chromosomes"));
       	addDialogComponent(new DialogComponentString(chromosome, ""));
    	setHorizontalPlacement(false);
    	/**General change listener**/
    	method.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(method.getStringValue().equals("geneanno")){
					exonsort.setEnabled(true);
					transcript_function.setEnabled(true);
					usesplicingthreshold.setEnabled(true);
					useindelsplicingthreshold.setEnabled(true);
					precedence.setEnabled(true);

					usescorethreshold.setBooleanValue(false);
					usescorethreshold.setEnabled(false);
		
				}
				else{
					exonsort.setBooleanValue(false);
					exonsort.setEnabled(false);
					transcript_function.setEnabled(false);
					transcript_function.setBooleanValue(false);
					usescorethreshold.setEnabled(true);
					usesplicingthreshold.setEnabled(false);
					useindelsplicingthreshold.setEnabled(false);
					usesplicingthreshold.setBooleanValue(false);
					useindelsplicingthreshold.setBooleanValue(false);
					precedence.setEnabled(false);
					
					
				}
				if(method.getStringValue().equals("regionanno")){
					usegff3dbfile.setEnabled(true);
					usebedfile.setEnabled(true);
					usescorecolumn.setEnabled(true);
					usecolswanted.setEnabled(true);
					
					minqueryfrac.setEnabled(true);
					usenormscorethreshold.setEnabled(true);

				}
				else{
					usegff3dbfile.setEnabled(false);
					usegff3dbfile.setBooleanValue(false);
					usebedfile.setEnabled(false);
					usebedfile.setBooleanValue(false);
					usescorecolumn.setEnabled(false);
					usescorecolumn.setBooleanValue(false);
					usecolswanted.setEnabled(false);
					usecolswanted.setBooleanValue(false);
					minqueryfrac.setEnabled(false);
					usenormscorethreshold.setBooleanValue(false);
					usenormscorethreshold.setEnabled(false);
					
					
				}
				if(method.getStringValue().equals("filter")){
					usegenericdbfile.setEnabled(true);
					usevcfdbfile.setEnabled(true);
					otherinfo.setEnabled(true);
			    	infoasscore.setEnabled(true);
			    	
			    	maf_threshold.setEnabled(true);
			    	indexfilter_threshold.setEnabled(true);
			    	usesift_threshold.setEnabled(true);

				}
				else{
					usegenericdbfile.setEnabled(false);
					usegenericdbfile.setBooleanValue(false);
					usevcfdbfile.setBooleanValue(false);
					usevcfdbfile.setEnabled(false);
					otherinfo.setEnabled(false);
			    	infoasscore.setEnabled(false);
			    	otherinfo.setBooleanValue(false);
			    	infoasscore.setBooleanValue(false);
			    	
			    	maf_threshold.setEnabled(false);
			    	indexfilter_threshold.setEnabled(false);
			    	usesift_threshold.setEnabled(false);
			    	usesift_threshold.setBooleanValue(false);
				}
			}
		});
    	usetablename.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usetablename.getBooleanValue()){
					tablename.setEnabled(true);
					queryfile.setEnabled(false);
				}
				else{
					tablename.setEnabled(false);
					queryfile.setEnabled(true);
				}
			}
		});
    	usesift_threshold.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usesift_threshold.getBooleanValue()){
					sift_threshold.setEnabled(true);
				}
				else{
					sift_threshold.setEnabled(false);
				}
			}
		});
    	/**In/Out changeListeners**/
    	/*useoutfile.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(useoutfile.getBooleanValue()){
					outfile.setEnabled(true);
				}
				else{
					outfile.setEnabled(false);
				}
			}
		});*/
    	/*usebuildver.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usebuildver.getBooleanValue()){
					buildver.setEnabled(true);
				}
				else{
					buildver.setEnabled(false);
				}
			}
		});
    	usedbtype.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usedbtype.getBooleanValue()){
					dbtype.setEnabled(true);
				}
				else{
					dbtype.setEnabled(false);
				}
			}
		});*/
    	usegff3dbfile.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usegff3dbfile.getBooleanValue()){
					gff3dbfile.setEnabled(true);
				}
				else{
					gff3dbfile.setEnabled(false);
				}
			}
		});
    	usebedfile.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usebedfile.getBooleanValue()){
					bedfile.setEnabled(true);
				}
				else{
					bedfile.setEnabled(false);
				}
			}
		});
    	usecolswanted.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usecolswanted.getBooleanValue()){
					colswanted.setEnabled(true);
				}
				else{
					colswanted.setEnabled(false);
				}
			}
		});
    	usescorecolumn.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usescorecolumn.getBooleanValue()){
					scorecolumn.setEnabled(true);
				}
				else{
					scorecolumn.setEnabled(false);
				}
			}
		});
    	usegenericdbfile.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usegenericdbfile.getBooleanValue()){
					genericdbfile.setEnabled(true);
				}
				else{
					genericdbfile.setEnabled(false);
				}
			}
		});
    	usevcfdbfile.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usevcfdbfile.getBooleanValue()){
					vcfdbfile.setEnabled(true);
				}
				else{
					vcfdbfile.setEnabled(false);
				}
			}
		});
		
    	
    	/**Finetuning changeListeners**/
    	useneargene.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(useneargene.getBooleanValue()){
					neargene.setEnabled(true);
					reverse.setEnabled(true);
				}
				else{
					neargene.setEnabled(false);
					reverse.setEnabled(false);
					reverse.setBooleanValue(false);
				}
			}
		});
    	usenormscorethreshold.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usenormscorethreshold.getBooleanValue()){
					normscore_threshold.setEnabled(true);
				}
				else{
					normscore_threshold.setEnabled(false);
				}
			}
		});
    	usescorethreshold.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usescorethreshold.getBooleanValue()){
					score_threshold.setEnabled(true);

				}
				else{
					score_threshold.setEnabled(false);
				}
			}
		});
    	usesplicingthreshold.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usesplicingthreshold.getBooleanValue()){
					splicing_threshold.setEnabled(true);
					indel_splicing_threshold.setEnabled(false);
					useindelsplicingthreshold.setBooleanValue(false);
				}
				else{
					splicing_threshold.setEnabled(false);
				}
			}
		});
    	useindelsplicingthreshold.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(useindelsplicingthreshold.getBooleanValue()){
					indel_splicing_threshold.setEnabled(true);
					usesplicingthreshold.setBooleanValue(false);
				}
				else{
					indel_splicing_threshold.setEnabled(false);
				}
			}
		});
    	usechromosome.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usechromosome.getBooleanValue()){
					chromosome.setEnabled(true);
				}
				else{
					chromosome.setEnabled(false);
				}
			}
		});

    }
}

