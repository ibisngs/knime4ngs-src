package de.helmholtz_muenchen.ibis.ngs.getannotationdatabase;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "GetAnnotationDatabase" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 */
public class GetAnnotationDatabaseNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring the GetAnnotationDatabase node.
     */
    protected GetAnnotationDatabaseNodeDialog() {
    	
    	final SettingsModelString organism = new SettingsModelString(GetAnnotationDatabaseNodeModel.CFGKEY_ORGANISM,"");
    	final SettingsModelString database = new SettingsModelString(GetAnnotationDatabaseNodeModel.CFGKEY_DATABASE,"");
    	final SettingsModelString annoType = new SettingsModelString(GetAnnotationDatabaseNodeModel.CFGKEY_ANNOTATIONTYPE,"");
    	final SettingsModelString outname = new SettingsModelString(GetAnnotationDatabaseNodeModel.CFGKEY_OUTNAME,"");

    	database.setEnabled(false);
    	outname.setEnabled(false);
    	outname.setStringValue("humandb");
    	database.setStringValue("hg18");
    	
    	createNewGroup("Annovar installpath");
    	addDialogComponent(new DialogComponentFileChooser(new SettingsModelString(GetAnnotationDatabaseNodeModel.CFGKEY_INSTALLPATH,null), "his_annodb_id", 0, true));
    	createNewGroup("General Options");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentStringSelection(organism,"Select organism", "Human", "Chimpanzee", "Macaque", "Cow", "Yeast", "Other DB"));
    	addDialogComponent(new DialogComponentString(database, "Name of database:", true, 7));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentString(outname, "Name of storage folder (any name):", true, 7));
    	addDialogComponent(new DialogComponentStringSelection(annoType,"Select annotation type", "[UCSC:gene-based] RefSeq transcript annotations and mRNA sequences in FASTA formats", 
    			"[UCSC:gene-based] UCSC Gene annotation (more comprehensive than RefSeq annotation) and mRNA sequences in FASTA format", 
    			"[UCSC:gene-based] Ensembl Gene annotation (more comprehensive than RefSeq annotation) and mRNA sequences in FASTA format", 
    			"[UCSC:gene-based] GENCODE manual annotation for hg18 (human)", "[UCSC:gene-based] GENCODE manual annotation for hg19 (human)",
    			"[UCSC:region-based] Approximate location of bands seen on Giemsa-stained chromosomes", 
    			"[UCSC:region-based] Transcription factor binding sites conserved in the human/mouse/rat alignment, based on transfac Matrix Database", 
    			"[UCSC:region-based] snoRNA and miRNA annotations", 
    			"[UCSC:region-based] TargetScan generated miRNA target site predictions", 
    			"[UCSC:region-based] Segmental duplications in genome", 
    			"[UCSC:region-based] Conserved elements produced by the phastCons program based on a whole-genome alignment of vertebrates", 
    			"[UCSC:region-based] Conserved functional RNA, through RNA secondary structure predictions made with the EvoFold program", 
    			"[UCSC:region-based] Database of Genomic Variants, which contains annotations for reported structural variations", 
    			"[UCSC:region-based] Canonical UCSC genes that have been associated with identifiers in the OMIM database", 
    			"[UCSC:region-based] Published GWAS results on diverse human diseases", 
    			"[1000 Genomes Project] Pilot 1 allele frequency data on the CEU, YRI and JPTCHB populations, updated on April 2009", 
    			"[1000 Genomes Project] Pilot 1 allele frequency data on the CEU, YRI and JPTCHB populations, updated on March 2010", 
    			"[1000 Genomes Project] Pilot 1 allele frequency data on the CEU, YRI and JPTCHB populations, updated on July 2010", 
    			"[1000 Genomes Project] Full phase 1 project with 629 subjects form diverse populations using August 2010, updated on November 2010", 
    			"[1000 Genomes Project] PHASE 1 low-coverage data on 1094 subjects using November 2010 alignments, updated on May 2011", 
    			"[1000 Genomes Project] variant calls from 1092 samples for SNPs, short indels, updated on February 2012", 
    			"[dbSNP] Version 132", 
    			"[dbSNP] Version 131", 
    			"[dbSNP] Version 130", 
    			"[dbSNP] Version 129", 
    			"[dbSNP] version 125", 
    			"[SIFT] Chr, start, end, reference allele, observed allele, SIFT score, reference amino acid, observed amino acids, updated on March 2011", 
    			"[dbNSFP] Lightweight database of human nonsynonymous SNPs and their functional predictions (PolyPhen2)", 
    			"[dbNSFP] Lightweight database of human nonsynonymous SNPs and their functional predictions (using LJBSIFT scores)", 
    			"[dbNSFP] Lightweight database of human nonsynonymous SNPs and their functional predictions (MutationTaster)", 
    			"[dbNSFP] Lightweight database of human nonsynonymous SNPs and their functional predictions (PhyloP conservation score)", 
    			"[dbNSFP] Lightweight database of human nonsynonymous SNPs and their functional predictions (LRT)", 
    			"[dbNSFP] Lightweight database of human nonsynonymous SNPs and their functional predictions (GERP++ score for exonic variants)", 
    			"[dbNSFP] Lightweight database of human nonsynonymous SNPs and their functional predictions (all scores above from LJB database)", 
    			"[ESP] 5400 NHLBI exomes (European Americans)", 
    			"[ESP] 5400 NHLBI exomes (African Americans)", 
    			"[ESP] 5400 NHLBI exomes (all ethnicity)"));
    	
    	
    	organism.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				database.setEnabled(organism.getStringValue().equals("Other DB"));
				outname.setEnabled(organism.getStringValue().equals("Other DB"));
				if(organism.getStringValue().equals("Human")) {
					outname.setStringValue("humandb");
					database.setStringValue("hg18");
				} else if(organism.getStringValue().equals("Chimpanzee")) {
					outname.setStringValue("chimpdb");
					database.setStringValue("panTro2");
				} else if(organism.getStringValue().equals("Macaque")) {
					outname.setStringValue("macaquedb");
					database.setStringValue("rheMac2");
				} else if(organism.getStringValue().equals("Yeast")) {
					outname.setStringValue("yeastdb");
					database.setStringValue("sacCer2");
				} else if(organism.getStringValue().equals("Cow")) {
					outname.setStringValue("cowdb");
					database.setStringValue("bosTau6");
				}
			}
		});
    	
    	annoType.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				String t = annoType.getStringValue();
				if(t.equals("[1000 Genomes Project] Pilot 1 allele frequency data on the CEU, YRI and JPTCHB populations, updated on April 2009") ||
						t.equals("[1000 Genomes Project] Pilot 1 allele frequency data on the CEU, YRI and JPTCHB populations, updated on March 2010") ||
						t.equals("[1000 Genomes Project] Pilot 1 allele frequency data on the CEU, YRI and JPTCHB populations, updated on July 2010") ||
						t.equals("[1000 Genomes Project] Full phase 1 project with 629 subjects form diverse populations using August 2010, updated on November 2010") ||
						t.equals("[1000 Genomes Project] PHASE 1 low-coverage data on 1094 subjects using November 2010 alignments, updated on May 2011") ||
						t.equals("[1000 Genomes Project] variant calls from 1092 samples for SNPs, short indels, updated on February 2012") ||
						t.equals("[dbSNP] Version 132") ||
						t.equals("[dbSNP] Version 131") ||
						t.equals("[dbSNP] Version 130") ||
						t.equals("[dbSNP] Version 129") ||
						t.equals("[dbSNP] version 125") ||
						t.equals("[SIFT] Chr, start, end, reference allele, observed allele, SIFT score, reference amino acid, observed amino acids, updated on March 2011") ||
						t.equals("[dbNSFP] Lightweight database of human nonsynonymous SNPs and their functional predictions (PolyPhen2)") ||
						t.equals("[dbNSFP] Lightweight database of human nonsynonymous SNPs and their functional predictions (using LJBSIFT scores)") ||
						t.equals("[dbNSFP] Lightweight database of human nonsynonymous SNPs and their functional predictions (MutationTaster)") ||
						t.equals("[dbNSFP] Lightweight database of human nonsynonymous SNPs and their functional predictions (PhyloP conservation score)") ||
						t.equals("[dbNSFP] Lightweight database of human nonsynonymous SNPs and their functional predictions (LRT)") ||
						t.equals("[dbNSFP] Lightweight database of human nonsynonymous SNPs and their functional predictions (GERP++ score for exonic variants)") ||
						t.equals("[dbNSFP] Lightweight database of human nonsynonymous SNPs and their functional predictions (all scores above from LJB database)") ||
						t.equals("[ESP] 5400 NHLBI exomes (European Americans)") ||
						t.equals("[ESP] 5400 NHLBI exomes (African Americans)") ||
						t.equals("[ESP] 5400 NHLBI exomes (all ethnicity)")) {
					organism.setStringValue("Human");
					outname.setStringValue("humandb");
			    	database.setStringValue("hg18");
			    	database.setEnabled(false);
			    	outname.setEnabled(false);
				}
			}
		});
    	
    }
}

