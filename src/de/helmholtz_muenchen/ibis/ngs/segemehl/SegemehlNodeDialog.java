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
package de.helmholtz_muenchen.ibis.ngs.segemehl;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "Segemehl" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Jan Quell
 * @author Maximilian Hastreiter
 * 
 */
public class SegemehlNodeDialog extends HTExecutorNodeDialog {

	
	
    /**
     * New pane for configuring the Segemehl node.
     */
    protected SegemehlNodeDialog() {}
    
    public void addToolDialogComponents() {
    	
    	final SettingsModelString segemehlfile = new SettingsModelString(SegemehlNodeModel.CFGKEY_SEGEMEHLFILE,"");
    	final SettingsModelString refseq = new SettingsModelString(SegemehlNodeModel.CFGKEY_REFSEQFILE,null);
    	
    	final SettingsModelBoolean clip5adapter = new SettingsModelBoolean(SegemehlNodeModel.CFGKEY_CLIP5ADAPTER, false);
    	final SettingsModelBoolean clip3adapter = new SettingsModelBoolean(SegemehlNodeModel.CFGKEY_CLIP3ADAPTER, false);
    	final SettingsModelBoolean autoadapter3seq = new SettingsModelBoolean(SegemehlNodeModel.CFGKEY_AUTOADAPTER3SEQ, false);
    	final SettingsModelString adapter3seq = new SettingsModelString(SegemehlNodeModel.CFGKEY_ADAPTER3SEQ, "");
    	final SettingsModelString adapter5seq = new SettingsModelString(SegemehlNodeModel.CFGKEY_ADAPTER5SEQ, "");
    	final SettingsModelBoolean clippolya = new SettingsModelBoolean(SegemehlNodeModel.CFGKEY_CLIPPOLYA, false);
    	final SettingsModelIntegerBounded clippingaccuracy = new SettingsModelIntegerBounded(SegemehlNodeModel.CFGKEY_CLIPPINGACCURACY, 70, 0, 100);
    	final SettingsModelString softhardclipping= new SettingsModelString(SegemehlNodeModel.CFGKEY_SOFTHARDCLIPPING,"Soft (Default)");
    	final SettingsModelBoolean checkBisulfiteMapping = new SettingsModelBoolean(SegemehlNodeModel.CFGKEY_CHECKSBISULFITEMAPPING, false);
    	final SettingsModelString bisulfiteMappingType = new SettingsModelString(SegemehlNodeModel.CFGKEY_BISULFITEMAPPINGTYPE,"methylC-seq/Lister et al.");
    	final SettingsModelIntegerBounded accuracy = new SettingsModelIntegerBounded(SegemehlNodeModel.CFGKEY_ACCURACY, 90, 0, 100);
    	final SettingsModelOptionalString optional = new SettingsModelOptionalString(SegemehlNodeModel.CFGKEY_OPTIONAL, "",false);
    	
    	addPrefPageSetting(segemehlfile, IBISKNIMENodesPlugin.SEGEMEHL);
    	addPrefPageSetting(refseq, IBISKNIMENodesPlugin.REF_GENOME);
    	
    	autoadapter3seq.setEnabled(false);
    	adapter3seq.setEnabled(false);
    	adapter5seq.setEnabled(false);
    	clippingaccuracy.setEnabled(false);
    	softhardclipping.setEnabled(false);
    	bisulfiteMappingType.setEnabled(false);
//    	checkBisulfiteMapping.setEnabled(false);
    	
    	createNewGroup("Alignment parameters");
    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(SegemehlNodeModel.CFGKEY_CHECKSPLITREADMAPPING, false), "Use (multiple) split read mapping (e.g. for cDNA reads)"));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(checkBisulfiteMapping, "Use bisulfite mapping"));
    	addDialogComponent(new DialogComponentStringSelection(bisulfiteMappingType,"with:","methylC-seq/Lister et al.","bs-seq/Cokus et al. protocol","PAR-CLIP with 4SU","PAR-CLIP with 6SG"));
    	addDialogComponent(new DialogComponentNumber(new SettingsModelIntegerBounded(SegemehlNodeModel.CFGKEY_THREADS, 4, 1, 250), "Number of threads/ cores to use:", 1));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentNumber(accuracy, "Alignment accuracy [%]:", 1));
    	addDialogComponent(new DialogComponentOptionalString(optional, "Optional Mapping Parameters"));
    	createNewGroup("Adapter and polyA clipping");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(clip5adapter, "Clip 5' adapters"));
    	addDialogComponent(new DialogComponentString(adapter5seq, "Adapter sequence:", true, 10));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(clip3adapter, "Clip 3' adapters"));
    	addDialogComponent(new DialogComponentString(adapter3seq, "Adapter sequence:", true, 10));
//    	addDialogComponent(new DialogComponentBoolean(autoadapter3seq, "Automatic 3’ adapter detection"));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentBoolean(clippolya, "Automatically clip polyA tails"));
    	addDialogComponent(new DialogComponentNumber(clippingaccuracy, "Clipping accuracy [%]:", 1));
    	addDialogComponent(new DialogComponentStringSelection(softhardclipping,"Type of clipping:","Soft (Default)","Hard"));
    	

    	
    	clip5adapter.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				adapter5seq.setEnabled(clip5adapter.getBooleanValue());
				clippingaccuracy.setEnabled(clippolya.getBooleanValue() || clip5adapter.getBooleanValue() || clip3adapter.getBooleanValue());
				softhardclipping.setEnabled(clippolya.getBooleanValue() || clip5adapter.getBooleanValue() || clip3adapter.getBooleanValue());
			}
		});

    	clip3adapter.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				adapter3seq.setEnabled(clip3adapter.getBooleanValue());
				autoadapter3seq.setEnabled(clip3adapter.getBooleanValue());
				clippingaccuracy.setEnabled(clippolya.getBooleanValue() || clip5adapter.getBooleanValue() || clip3adapter.getBooleanValue());
				softhardclipping.setEnabled(clippolya.getBooleanValue() || clip5adapter.getBooleanValue() || clip3adapter.getBooleanValue());
			}
		});
    	
    	autoadapter3seq.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				adapter3seq.setEnabled(clip3adapter.getBooleanValue() && !autoadapter3seq.getBooleanValue());
			}
		});
    	
    	clippolya.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				clippingaccuracy.setEnabled(clippolya.getBooleanValue() || clip5adapter.getBooleanValue() || clip3adapter.getBooleanValue());
				softhardclipping.setEnabled(clippolya.getBooleanValue() || clip5adapter.getBooleanValue() || clip3adapter.getBooleanValue());
			}
		});
    	
    	checkBisulfiteMapping.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				bisulfiteMappingType.setEnabled(checkBisulfiteMapping.getBooleanValue());
			}
		});

    }	
}

