package de.helmholtz_muenchen.ibis.ngs.depthofcoverage;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeDialog;

/**
 * <code>NodeDialog</code> for the "DepthOfCoverage" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class DepthOfCoverageNodeDialog extends GATKNodeDialog {

	
//	private final SettingsModelOptionalString m_path2bed = new SettingsModelOptionalString(DepthOfCoverageNodeModel.CFGKEY_PATH2BED,"",true);
//	private final SettingsModelBoolean m_bed_file_check = new SettingsModelBoolean(DepthOfCoverageNodeModel.CFGKEY_BED_FILE_CHECKBOX,false);
//	private final SettingsModelOptionalString m_extrafilters = new SettingsModelOptionalString(DepthOfCoverageNodeModel.CFGKEY_EXTRAFILTERS,"",false);
//	private final SettingsModelIntegerBounded m_threads = new SettingsModelIntegerBounded(DepthOfCoverageNodeModel.CFGKEY_THREADS, 1, 1,Integer.MAX_VALUE);
//	private final SettingsModelString m_filesuffix = new SettingsModelString(DepthOfCoverageNodeModel.CFGKEY_FILESUFFIX,"");
//	private final SettingsModelString m_infile = new SettingsModelString(DepthOfCoverageNodeModel.CFGKEY_INFILE,"");
	
	
    /**
     * New pane for configuring the DepthOfCoverage node.
     */
//    protected DepthOfCoverageNodeDialog() {
//
//    	createNewGroup("Path to GATK jar file");
//    	addDialogComponent(new DialogComponentFileChooser(m_path2gatk, "his_id_GATK_DoC", 0, ".jar"));
//    	createNewGroup("Path to Infile");
//    	addDialogComponent(new DialogComponentFileChooser(m_infile, "his_id_GATK_DoC", 0, ".bam"));
//    	createNewGroup("Path to Reference Genome");
//    	addDialogComponent(new DialogComponentFileChooser(m_refgenome, "his_id_GATK_DoC", 0, ".fasta|.fa"));
//    	createNewGroup("Path to BED file");
//    	addDialogComponent(new DialogComponentFileChooser(m_path2bed, "his_id_GATK_DoC", 0, ".bed"));
//    	
//    	addDialogComponent(new DialogComponentString(m_filesuffix, "Outfile Suffix"));
//    	addDialogComponent(new DialogComponentOptionalString(m_extrafilters, "GATK Filter Terms (separated by ',')"));
//    	
//    	addDialogComponent(new DialogComponentNumber(m_memory, "GATK Memory", 1));
//    	addDialogComponent(new DialogComponentNumber(m_threads, "GATK Threads", 1));
//    	
//    	
//    }


	@Override
	protected void addDialogComponent() {

		final SettingsModelString m_path2bed = new SettingsModelOptionalString(DepthOfCoverageNodeModel.CFGKEY_PATH2BED,"",true);
		final SettingsModelBoolean m_bed_file_check = new SettingsModelBoolean(DepthOfCoverageNodeModel.CFGKEY_BED_FILE_CHECKBOX,false);
		final SettingsModelOptionalString m_extrafilters = new SettingsModelOptionalString(DepthOfCoverageNodeModel.CFGKEY_EXTRAFILTERS,"",false);
		final SettingsModelString m_filesuffix = new SettingsModelString(DepthOfCoverageNodeModel.CFGKEY_FILESUFFIX,"");
		final SettingsModelString m_infile = new SettingsModelString(DepthOfCoverageNodeModel.CFGKEY_INFILE,"");
		
		createNewTab("DoC options");
    	createNewGroup("Path to Infile");
    	addDialogComponent(new DialogComponentFileChooser(m_infile, "his_id_GATK_DoC", 0, ".bam"));

    	createNewGroup("Path to BED file");
    	addDialogComponent(new DialogComponentBoolean(m_bed_file_check,"Use bed file?"));
    	addDialogComponent(new DialogComponentFileChooser(m_path2bed, "his_id_GATK_DoC", 0, ".bed"));
    	
    	addDialogComponent(new DialogComponentString(m_filesuffix, "Outfile Suffix"));
    	addDialogComponent(new DialogComponentOptionalString(m_extrafilters, "GATK Filter Terms (separated by ',')"));
		
	}
}

