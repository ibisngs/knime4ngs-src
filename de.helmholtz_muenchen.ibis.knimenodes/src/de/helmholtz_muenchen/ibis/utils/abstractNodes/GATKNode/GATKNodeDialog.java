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
package de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.ngs.depthofcoverage.DepthOfCoverageNodeModel;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;



public abstract class GATKNodeDialog extends HTExecutorNodeDialog{

//    private final SettingsModelString GATK = new SettingsModelString(GATKNodeModel.CFGKEY_GATK_PATH, "");
//    private final SettingsModelString REF_GENOME = new SettingsModelString(GATKNodeModel.CFGKEY_REF_GENOME, "");
//    private final SettingsModelIntegerBounded m_GATK_MEM = new SettingsModelIntegerBounded(GATKNodeModel.CFGKEY_GATK_MEM, 4, 1, Integer.MAX_VALUE);
//	  private final SettingsModelString m_path2bed = new SettingsModelOptionalString(DepthOfCoverageNodeModel.CFGKEY_PATH2BED,"",true);
//	  private final SettingsModelBoolean m_bed_file_check = new SettingsModelBoolean(DepthOfCoverageNodeModel.CFGKEY_BED_FILE_CHECKBOX,false);
//    private final SettingsModelOptionalString m_OPT_FLAGS = new SettingsModelOptionalString(GATKNodeModel.CFGKEY_OPT_FLAGS,"",false);

    protected GATKNodeDialog() {
    }
    
    public void addToolDialogComponents() {
    	addDialogComponent();
    	
    	if(optionsTabIsEmpty()) {
			setDefaultTabTitle("GATK");
		} else {
			createNewTab("GATK");
		}
    	
    	final SettingsModelString GATK = new SettingsModelString(GATKNodeModel.CFGKEY_GATK_PATH, "");
        final SettingsModelString REF_GENOME = new SettingsModelString(GATKNodeModel.CFGKEY_REF_GENOME, "");
        final SettingsModelIntegerBounded m_GATK_MEM = new SettingsModelIntegerBounded(GATKNodeModel.CFGKEY_GATK_MEM, 4, 1, Integer.MAX_VALUE);
    	final SettingsModelString m_path2bed = new SettingsModelOptionalString(DepthOfCoverageNodeModel.CFGKEY_PATH2BED,"",true);
    	final SettingsModelBoolean m_bed_file_check = new SettingsModelBoolean(DepthOfCoverageNodeModel.CFGKEY_BED_FILE_CHECKBOX,false);
        final SettingsModelOptionalString m_OPT_FLAGS = new SettingsModelOptionalString(GATKNodeModel.CFGKEY_OPT_FLAGS,"",false);
    	
    	addPrefPageSetting(REF_GENOME, IBISKNIMENodesPlugin.REF_GENOME);
    	addPrefPageSetting(GATK, IBISKNIMENodesPlugin.GATK);
    	
    	addDialogComponent(new DialogComponentNumber(m_GATK_MEM, "GATK Memory", 1));
    	addDialogComponent(new DialogComponentBoolean(m_bed_file_check,"Use BED file?"));
    	m_path2bed.setEnabled(false);
    	
    	createNewGroup("Path to BED file");
    	addDialogComponent(new DialogComponentFileChooser(m_path2bed, "his_id_GATK_DoC", 0, ".bed"));
    	
    	createNewGroup("Further options");
    	addDialogComponent(new DialogComponentOptionalString(m_OPT_FLAGS,"Optional flags"));
    	
    	
    	m_bed_file_check.addChangeListener(new ChangeListener () {
    		
    		@Override
			public void stateChanged(ChangeEvent e) {
				m_path2bed.setEnabled(m_bed_file_check.getBooleanValue());
			}
    	});

    }
	
    protected abstract void addDialogComponent();
    
}
