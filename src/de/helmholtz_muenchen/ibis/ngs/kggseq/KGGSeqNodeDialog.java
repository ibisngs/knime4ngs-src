package de.helmholtz_muenchen.ibis.ngs.kggseq;

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
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;



/**
 * <code>NodeDialog</code> for the "KGGSeq" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class KGGSeqNodeDialog extends DefaultNodeSettingsPane {

	
    private final SettingsModelString m_KGGSEQ = new SettingsModelString(KGGSeqNodeModel.CFGKEY_KGGSEQ_PATH, "");
    private final SettingsModelString m_BUILDVER = new SettingsModelString(KGGSeqNodeModel.CFGKEY_BUILDVER, "hg19");
    private final SettingsModelString m_INFILE = new SettingsModelString(KGGSeqNodeModel.CFGKEY_INFILE, "");
    private final SettingsModelString m_PEDFILE = new SettingsModelString(KGGSeqNodeModel.CFGKEY_PEDFILE, "");
    private final SettingsModelBoolean m_COMPOSITESUBJECTID = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_COMPOSITE_SUBJECT_ID, false);
    private final SettingsModelString m_OUTPREFIX = new SettingsModelString(KGGSeqNodeModel.CFGKEY_OUTPREFIX, "");
    private final SettingsModelString m_OUTFORMAT = new SettingsModelString(KGGSeqNodeModel.CFGKEY_OUTFORMAT, "excel");
   
    private final SettingsModelDoubleBounded m_SEQ_QUAL = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_SEQ_QUAL,50,0,Double.MAX_VALUE);
    private final SettingsModelDoubleBounded m_SEQ_MQ = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_SEQ_MQ,20,0,Double.MAX_VALUE);
    private final SettingsModelDoubleBounded m_SEQ_SB = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_SEQ_SB,-10,-Double.MAX_VALUE,Double.MAX_VALUE);
    private final SettingsModelDoubleBounded m_GTY_QUAL = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_GTY_QUAL,20,0,Double.MAX_VALUE);
    private final SettingsModelDoubleBounded m_GTY_DP = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_GTY_DP,8,0,Double.MAX_VALUE);
    private final SettingsModelIntegerBounded m_GTY_SEC_PL = new SettingsModelIntegerBounded(KGGSeqNodeModel.CFGKEY_GTY_SEC_PL,20,0,Integer.MAX_VALUE);
    private final SettingsModelDoubleBounded m_GTY_AF_REF = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_GTY_AF_REF,0.05,0,Double.MAX_VALUE);
    private final SettingsModelDoubleBounded m_GTY_AF_HET = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_GTY_AF_HET,0.25,0,Double.MAX_VALUE);
    private final SettingsModelDoubleBounded m_GTY_AF_ALT = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_GTY_AF_ALT,0.75,0,Double.MAX_VALUE);

    private final SettingsModelOptionalString m_GENOTYPE_FILTER = new SettingsModelOptionalString(KGGSeqNodeModel.CFGKEY_GENOTYPE_FILTER, "4",true);
    private final SettingsModelBoolean m_IGNORE_HOMO = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_IGNORE_HOMO, true);
    private final SettingsModelOptionalString m_GENE_FEATURES = new SettingsModelOptionalString(KGGSeqNodeModel.CFGKEY_GENE_FEATURES, "0,1,2,3,4,5",true);
    private final SettingsModelBoolean m_FILTER_COMMON = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_FILTER_COMMON, true);   
    private final SettingsModelBoolean m_DISEASE_CAUSING_PRED = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_DISEASE_CAUSING_PRED, true);
    private final SettingsModelBoolean m_OMIM_ANNO = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_OMIM_ANNO, true);
    private final SettingsModelBoolean m_CANDIDATE_PPI = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_CANDIDATE_PPI, true);
    private final SettingsModelString m_CANDIDATE_GENES = new SettingsModelString(KGGSeqNodeModel.CFGKEY_CANDIDATE_GENES, "");
    private final SettingsModelBoolean m_CANDIDATE_PATHWAYS = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_CANDIDATE_PATHWAYS, true);
    private final SettingsModelOptionalString m_PUBMED = new SettingsModelOptionalString(KGGSeqNodeModel.CFGKEY_PUBMED,"",false);
	
	
	
    /**
     * New pane for configuring the KGGSeq node.
     */
    protected KGGSeqNodeDialog() {

    	createNewGroup("KGGSeq Jar");
    	addDialogComponent(new DialogComponentFileChooser(m_KGGSEQ, "kggseq", 0, ".jar"));
    	
    	createNewGroup("VCF Input File");
    	addDialogComponent(new DialogComponentFileChooser(m_INFILE, "kggseq_vcfin", 0, ".vcf")); 
    	
    	createNewGroup("PED Input File");
    	addDialogComponent(new DialogComponentFileChooser(m_PEDFILE, "kggseq_pedin", 0, ".ped")); 
    	
    	createNewGroup("General Options");
    	addDialogComponent(new DialogComponentStringSelection(m_BUILDVER, "Build Version", "hg19"));
    	addDialogComponent(new DialogComponentBoolean(m_COMPOSITESUBJECTID, "Composite Subject IDs"));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentString(m_OUTPREFIX, "Outfile Prefix"));
    	addDialogComponent(new DialogComponentStringSelection(m_OUTFORMAT, "Output Format","excel"));
    	setHorizontalPlacement(false);
    	
    	createNewGroup("Quality Cutoffs");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentNumber(m_SEQ_QUAL, "Seq Qual", 1,6));
    	addDialogComponent(new DialogComponentNumber(m_SEQ_MQ, "Seq MQ", 1,6));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentNumber(m_SEQ_SB, "Seq SB", 1,6));
    	addDialogComponent(new DialogComponentNumber(m_GTY_QUAL, "Gty Qual", 1,6));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentNumber(m_GTY_DP, "Gty DP", 1,6));
    	addDialogComponent(new DialogComponentNumber(m_GTY_SEC_PL, "Gty Sec Pl", 1,6));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentNumber(m_GTY_AF_REF, "Gty Af Ref", 0.05,6));
    	addDialogComponent(new DialogComponentNumber(m_GTY_AF_HET, "Gty Af Het", 0.05,6));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentNumber(m_GTY_AF_ALT, "Gty Af Alt", 0.05,6));
    	
    	createNewGroup("Annotation Options");
    	addDialogComponent(new DialogComponentOptionalString(m_GENOTYPE_FILTER, "Genotype Filter"));
    	addDialogComponent(new DialogComponentBoolean(m_IGNORE_HOMO, "Ignore Homo"));
    	addDialogComponent(new DialogComponentOptionalString(m_GENE_FEATURES, "Select Gene Features"));
    	addDialogComponent(new DialogComponentBoolean(m_FILTER_COMMON, "Filter by Common variants")); 
    	addDialogComponent(new DialogComponentBoolean(m_DISEASE_CAUSING_PRED, "Prioritize sequence variants by disease-causing prediction")); 
    	addDialogComponent(new DialogComponentBoolean(m_OMIM_ANNO, "Prioritize sequence variants by other genomic and OMIM annotation")); 
    	addDialogComponent(new DialogComponentOptionalString(m_PUBMED, "Prioritize sequence variants by PubMed"));
    	addDialogComponent(new DialogComponentBoolean(m_CANDIDATE_PPI, "Prioritize sequence variants by candidate genes with  protein interaction information")); 
    	addDialogComponent(new DialogComponentBoolean(m_CANDIDATE_PATHWAYS, "Prioritize sequence variants by candidate genes with pathway information")); 

    	createNewGroup("Candidate Genes File");
    	addDialogComponent(new DialogComponentFileChooser(m_CANDIDATE_GENES, "kggseq_candidate_genes", 0, ""));  
    	
    	
    	
    	m_CANDIDATE_PPI.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
					if(m_CANDIDATE_PPI.getBooleanValue() || m_CANDIDATE_PATHWAYS.getBooleanValue()){
						m_CANDIDATE_GENES.setEnabled(true);
					}else{
						m_CANDIDATE_GENES.setEnabled(false);
					}
			}
		});
    	m_CANDIDATE_PATHWAYS.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
					if(m_CANDIDATE_PPI.getBooleanValue() || m_CANDIDATE_PATHWAYS.getBooleanValue()){
						m_CANDIDATE_GENES.setEnabled(true);
					}else{
						m_CANDIDATE_GENES.setEnabled(false);
					}
			}
		});
    	
    }
}


