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

package de.helmholtz_muenchen.ibis.ngs.vepsummary;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;


/**
 * <code>NodeDialog</code> for the "LOFStatistics" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Tim Jeske
 */
public class VEPSummaryNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring the LOFStatistics node.
     */
	private final SettingsModelString cdsin = new SettingsModelString(VEPSummaryNodeModel.CFGKEY_CDS_INFILE,"-");
	private final SettingsModelString pedin = new SettingsModelString(VEPSummaryNodeModel.CFGKEY_PED_INFILE,"-");
	private final SettingsModelBoolean internal_gene_set = new SettingsModelBoolean(VEPSummaryNodeModel.CFGKEY_INTERNAL_GENE_SET,true);
	private final SettingsModelBoolean create_var_sum = new SettingsModelBoolean(VEPSummaryNodeModel.CFGKEY_CREATE_VAR_SUM,true);
//	private final SettingsModelBoolean create_gene_sum = new SettingsModelBoolean(LOFSummaryNodeModel.CFGKEY_CREATE_GENE_SUM,true);
	private final SettingsModelBoolean create_sample_sum = new SettingsModelBoolean(VEPSummaryNodeModel.CFGKEY_CREATE_SAMPLE_SUM,true);
	private final SettingsModelBoolean parallel_exec = new SettingsModelBoolean(VEPSummaryNodeModel.CFGKEY_PARALLEL_EXEC,false);
	private final SettingsModelBoolean create_matrix = new SettingsModelBoolean(VEPSummaryNodeModel.CFGKEY_CREATE_MATRIX,true);

	
    protected VEPSummaryNodeDialog() {
    	
    	createNewGroup("Path to CDS/GTF file");
    	addDialogComponent(new DialogComponentFileChooser(cdsin, "his_id_LOFStatistics_CDSIN", 0, ".fa|.fasta|.gtf"));
    	addDialogComponent(new DialogComponentBoolean(internal_gene_set, "Use internal representation of CDS/GTF file?"));
    	
    	createNewGroup("Path to PED file");
    	addDialogComponent(new DialogComponentFileChooser(pedin, "his_id_LOFStatistics_PEDIN", 0, ".ped"));
    	
    	
    	createNewGroup("Further options");
    	addDialogComponent(new DialogComponentBoolean(create_var_sum,"Create variant summary?"));
//    	addDialogComponent(new DialogComponentBoolean(create_gene_sum,"Create gene summary?"));
    	addDialogComponent(new DialogComponentBoolean(create_sample_sum,"Create sample summary?"));
    	addDialogComponent(new DialogComponentBoolean(create_matrix, "Create matrix file?"));
    	addDialogComponent(new DialogComponentBoolean(parallel_exec,"Create summaries in parallel?"));
    }
}

