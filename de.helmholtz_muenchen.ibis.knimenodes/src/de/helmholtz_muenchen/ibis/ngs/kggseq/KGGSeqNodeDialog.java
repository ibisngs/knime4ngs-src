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

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;



/**
 * <code>NodeDialog</code> for the "KGGSeq" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Maximilian Hastreiter
 */
public class KGGSeqNodeDialog extends HTExecutorNodeDialog {


    /**
     * New pane for configuring the KGGSeq node.
     */
    protected KGGSeqNodeDialog() {
    	super(VCFCell.TYPE.getPreferredValueClass(), 0);
    }

	@Override
	public void addToolDialogComponents() {
		
	    final SettingsModelString m_KGGSEQ = new SettingsModelString(KGGSeqNodeModel.CFGKEY_KGGSEQ_PATH, "");
	    final SettingsModelString m_BUILDVER = new SettingsModelString(KGGSeqNodeModel.CFGKEY_BUILDVER, "hg19");
	    final SettingsModelString m_PEDFILE = new SettingsModelString(KGGSeqNodeModel.CFGKEY_PEDFILE, "");
	    final SettingsModelString m_RESOURCE = new SettingsModelString(KGGSeqNodeModel.CFGKEY_RESOURCE, "");
	    final SettingsModelBoolean m_COMPOSITESUBJECTID = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_COMPOSITE_SUBJECT_ID, false);
	    final SettingsModelString m_OUTPREFIX = new SettingsModelString(KGGSeqNodeModel.CFGKEY_OUTPREFIX, "");
	    final SettingsModelString m_OUTFORMAT = new SettingsModelString(KGGSeqNodeModel.CFGKEY_OUTFORMAT, "excel");
	   
	    final SettingsModelDoubleBounded m_SEQ_QUAL = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_SEQ_QUAL,50,0,Double.MAX_VALUE);
	    final SettingsModelDoubleBounded m_SEQ_MQ = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_SEQ_MQ,20,0,Double.MAX_VALUE);
	    final SettingsModelDoubleBounded m_SEQ_SB = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_SEQ_SB,-10,-Double.MAX_VALUE,Double.MAX_VALUE);
	    final SettingsModelDoubleBounded m_GTY_QUAL = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_GTY_QUAL,20,0,Double.MAX_VALUE);
	    final SettingsModelDoubleBounded m_GTY_DP = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_GTY_DP,8,0,Double.MAX_VALUE);
	    final SettingsModelIntegerBounded m_GTY_SEC_PL = new SettingsModelIntegerBounded(KGGSeqNodeModel.CFGKEY_GTY_SEC_PL,20,0,Integer.MAX_VALUE);
	    final SettingsModelDoubleBounded m_GTY_AF_REF = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_GTY_AF_REF,0.05,0,Double.MAX_VALUE);
	    final SettingsModelDoubleBounded m_GTY_AF_HET = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_GTY_AF_HET,0.25,0,Double.MAX_VALUE);
	    final SettingsModelDoubleBounded m_GTY_AF_ALT = new SettingsModelDoubleBounded(KGGSeqNodeModel.CFGKEY_GTY_AF_ALT,0.75,0,Double.MAX_VALUE);

	    final SettingsModelOptionalString m_GENOTYPE_FILTER = new SettingsModelOptionalString(KGGSeqNodeModel.CFGKEY_GENOTYPE_FILTER, "4",true);
	    final SettingsModelBoolean m_IGNORE_HOMO = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_IGNORE_HOMO, false);
	    final SettingsModelOptionalString m_GENE_FEATURES = new SettingsModelOptionalString(KGGSeqNodeModel.CFGKEY_GENE_FEATURES, "0,1,2,3,4,5,6",true);
	    final SettingsModelOptionalString m_FILTER_COMMON = new SettingsModelOptionalString(KGGSeqNodeModel.CFGKEY_FILTER_COMMON,"hg19_1kg201204,hg19_dbsnp138,hg19_ESP6500AA,hg19_ESP6500EA", false);   
	    final SettingsModelBoolean m_DISEASE_CAUSING_PRED = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_DISEASE_CAUSING_PRED, false);
	    final SettingsModelBoolean m_OMIM_ANNO = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_OMIM_ANNO, false);
	    final SettingsModelBoolean m_CANDIDATE_PPI = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_CANDIDATE_PPI, false);
	    final SettingsModelString m_CANDIDATE_GENES = new SettingsModelString(KGGSeqNodeModel.CFGKEY_CANDIDATE_GENES, "");
	    final SettingsModelBoolean m_CANDIDATE_PATHWAYS = new SettingsModelBoolean(KGGSeqNodeModel.CFGKEY_CANDIDATE_PATHWAYS, false);
	    final SettingsModelOptionalString m_PUBMED = new SettingsModelOptionalString(KGGSeqNodeModel.CFGKEY_PUBMED,"",false);
		
	    addPrefPageSetting(m_KGGSEQ, IBISKNIMENodesPlugin.KGGSeq);
    	
    	createNewGroup("PED Input File");
    	addDialogComponent(new DialogComponentFileChooser(m_PEDFILE, "kggseq_pedin", 0, ".ped")); 
    	
    	createNewGroup("Local path to the resource datasets");
    	addDialogComponent(new DialogComponentFileChooser(m_RESOURCE, "kggseq_resource", 0,true)); 
    	
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
    	addDialogComponent(new DialogComponentOptionalString(m_FILTER_COMMON, "Filter by Common variants (AF=0.05)")); 
    	addDialogComponent(new DialogComponentBoolean(m_DISEASE_CAUSING_PRED, "Prioritize sequence variants by disease-causing prediction")); 
    	addDialogComponent(new DialogComponentBoolean(m_OMIM_ANNO, "Prioritize sequence variants by OMIM annotation")); 
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


