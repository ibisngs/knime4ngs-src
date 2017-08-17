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
package de.helmholtz_muenchen.ibis.ngs.gatkunifiedgenotyper;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeDialog;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.BAMCell;


/**
 * <code>NodeDialog</code> for the "GATKUnifiedGenotyper" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Marie-Sophie Friedl
 * @author Maximilian Hastreiter
 * @author Tim Jeske
 */
public class GATKUnifiedGenotyperNodeDialog extends GATKNodeDialog {
	
	public GATKUnifiedGenotyperNodeDialog() {
		super(BAMCell.TYPE.getPreferredValueClass());
	}
	
	public void addDialogComponent() {
		
		final SettingsModelString outfolder = new SettingsModelString(GATKUnifiedGenotyperNodeModel.CFGKEY_OUTFOLDER, "");
		final SettingsModelString variant_type= new SettingsModelString(GATKUnifiedGenotyperNodeModel.CFGKEY_VARIANT_TYPE, GATKUnifiedGenotyperNodeModel.DEF_VARIANT_TYPE);

		final SettingsModelBoolean use_dbsnp=new SettingsModelBoolean(GATKUnifiedGenotyperNodeModel.CFGKEY_USE_DBSNP, GATKUnifiedGenotyperNodeModel.DEF_USE_DBSNP);
		final SettingsModelString dbsnp_file=new SettingsModelString(GATKUnifiedGenotyperNodeModel.CFGKEY_DBSNP_FILE, GATKUnifiedGenotyperNodeModel.DEF_DBSNP_FILE);
		final SettingsModelIntegerBounded num_threads = new SettingsModelIntegerBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_NUM_THREADS, GATKUnifiedGenotyperNodeModel.DEF_NUM_THREADS, GATKUnifiedGenotyperNodeModel.MIN_NUM_THREADS, GATKUnifiedGenotyperNodeModel.MAX_NUM_THREADS);
	
		final SettingsModelBoolean mbq = new SettingsModelBoolean(GATKUnifiedGenotyperNodeModel.CFGKEY_MBQ, GATKUnifiedGenotyperNodeModel.DEF_MBQ);
		final SettingsModelDoubleBounded call_min_confidence= new SettingsModelDoubleBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_CALL_MIN_CONFIDENCE, GATKUnifiedGenotyperNodeModel.DEF_CALL_MIN_CONFIDENCE, GATKUnifiedGenotyperNodeModel.MIN_CALL_MIN_CONFIDENCE, GATKUnifiedGenotyperNodeModel.MAX_CALL_MIN_CONFIDENCE);
		final SettingsModelDoubleBounded emit_min_confidence = new SettingsModelDoubleBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_EMIT_MIN_CONFIDENCE, GATKUnifiedGenotyperNodeModel.DEF_EMIT_MIN_CONFIDENCE, GATKUnifiedGenotyperNodeModel.MIN_EMIT_MIN_CONFIDENCE, GATKUnifiedGenotyperNodeModel.MAX_EMIT_MIN_CONFIDENCE);
		final SettingsModelDoubleBounded pcr_error = new SettingsModelDoubleBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_PCR_ERR, GATKUnifiedGenotyperNodeModel.DEF_PCR_ERR, GATKUnifiedGenotyperNodeModel.MIN_PCR_ERR, GATKUnifiedGenotyperNodeModel.MAX_PCR_ERR);
		final SettingsModelDoubleBounded contamination = new SettingsModelDoubleBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_CONTAMINATION, GATKUnifiedGenotyperNodeModel.DEF_CONTAMINATION, GATKUnifiedGenotyperNodeModel.MIN_CONTAMINATION, GATKUnifiedGenotyperNodeModel.MAX_CONTAMINATION);
		final SettingsModelDoubleBounded heterozygosity= new SettingsModelDoubleBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_HET, GATKUnifiedGenotyperNodeModel.DEF_HET, GATKUnifiedGenotyperNodeModel.MIN_HET, GATKUnifiedGenotyperNodeModel.MAX_HET);
		final SettingsModelDoubleBounded max_deletion_fraction= new SettingsModelDoubleBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_MAX_DELETION_FRAC, GATKUnifiedGenotyperNodeModel.DEF_MAX_DELETION_FRAC, GATKUnifiedGenotyperNodeModel.MIN_MAX_DELETION_FRAC, GATKUnifiedGenotyperNodeModel.MAX_MAX_DELETION_FRAC);
		final SettingsModelIntegerBounded min_base_qual= new SettingsModelIntegerBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_MIN_BASE_QUAL, GATKUnifiedGenotyperNodeModel.DEF_MIN_BASE_QUAL, GATKUnifiedGenotyperNodeModel.MIN_MIN_BASE_QUAL, GATKUnifiedGenotyperNodeModel.MAX_MIN_BASE_QUAL);
		final SettingsModelDoubleBounded indel_heterozygosity= new SettingsModelDoubleBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_INDEL_HET, GATKUnifiedGenotyperNodeModel.DEF_INDEL_HET, GATKUnifiedGenotyperNodeModel.MIN_INDEL_HET, GATKUnifiedGenotyperNodeModel.MAX_INDEL_HET);
		final SettingsModelIntegerBounded min_indel_count = new SettingsModelIntegerBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_MIN_INDEL_CNT, GATKUnifiedGenotyperNodeModel.DEF_MIN_INDEL_CNT, GATKUnifiedGenotyperNodeModel.MIN_MIN_INDEL_CNT, GATKUnifiedGenotyperNodeModel.MAX_MIN_INDEL_CNT);
		final SettingsModelDoubleBounded min_indel_frac = new SettingsModelDoubleBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_MIN_INDEL_FRAC, GATKUnifiedGenotyperNodeModel.DEF_MIN_INDEL_FRAC, GATKUnifiedGenotyperNodeModel.MIN_MIN_INDEL_FRAC, GATKUnifiedGenotyperNodeModel.MAX_MIN_INDEL_FRAC);
		final SettingsModelIntegerBounded gap_open_pen= new SettingsModelIntegerBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_GAP_OPEN_PEN, GATKUnifiedGenotyperNodeModel.DEF_GAP_OPEN_PEN, GATKUnifiedGenotyperNodeModel.MIN_GAP_OPEN_PEN, GATKUnifiedGenotyperNodeModel.MAX_GAP_OPEN_PEN);
		final SettingsModelIntegerBounded gap_cont_pen= new SettingsModelIntegerBounded(GATKUnifiedGenotyperNodeModel.CFGKEY_GAP_CONT_PEN, GATKUnifiedGenotyperNodeModel.DEF_GAP_CONT_PEN, GATKUnifiedGenotyperNodeModel.MIN_GAP_CONT_PEN, GATKUnifiedGenotyperNodeModel.MAX_GAP_CONT_PEN);
	
		//Proxy options
		//	private final SettingsModelBoolean useproxy = new SettingsModelBoolean(GATKUnifiedGenotyperNodeModel.CFGKEY_USEPROXY, false);
		//	final SettingsModelString proxyhost = new SettingsModelString(GATKUnifiedGenotyperNodeModel.CFGKEY_PROXYHOST, null);
		//	final SettingsModelString proxyport = new SettingsModelString(GATKUnifiedGenotyperNodeModel.CFGKEY_PROXYPORT, null);
		//	private final SettingsModelBoolean useproxyauth = new SettingsModelBoolean(GATKUnifiedGenotyperNodeModel.CFGKEY_USEPROXYAUTH, false);
		//	final SettingsModelString proxyuser = new SettingsModelString(GATKUnifiedGenotyperNodeModel.CFGKEY_PROXYUSER, null);
		//	final SettingsModelString proxypassword = new SettingsModelString(GATKUnifiedGenotyperNodeModel.CFGKEY_PROXYPASSWORD, null);
	    
		addPrefPageSetting(dbsnp_file, IBISKNIMENodesPlugin.RES_DBSNP);

        addDialogComponent(new DialogComponentStringSelection(variant_type, "Choose variant type(s) to be called", GATKUnifiedGenotyperNodeModel.AVAIL_VARIANT_TYPE));
		addDialogComponent(new DialogComponentBoolean(use_dbsnp, "Use dbSNP file to populate ID column?"));
		
		createNewGroup("Output folder");
		addDialogComponent(new DialogComponentFileChooser(outfolder, "his_ID_outfolder",0, true));
         
        createNewGroup("SNP calling");
        addDialogComponent(new DialogComponentNumber(min_base_qual, "Minimum base quality score for calling", 1, 5));
        
        createNewGroup("Indel calling");
        addDialogComponent(new DialogComponentNumber(indel_heterozygosity, "Heterozygosity for indel calling",0.00001 , 5));
        addDialogComponent(new DialogComponentNumber(min_indel_count, "Minimum count of indel reads", 1, 5));
        addDialogComponent(new DialogComponentNumber(min_indel_frac, "Minimum fraction of indel reads ", 0.01, 5));
        addDialogComponent(new DialogComponentNumber(gap_open_pen, "Indel gap open penalty", 1, 5));
        addDialogComponent(new DialogComponentNumber(gap_cont_pen, "Indel continuation penalty", 1, 5));
        
        createNewGroup("Further parameters");
        addDialogComponent(new DialogComponentNumber(contamination, "Fraction of contamination in sequencing data", 0.0001, 5));
        addDialogComponent(new DialogComponentNumber(heterozygosity, "Heterozygosity", 0.0001, 5));
        addDialogComponent(new DialogComponentNumber(max_deletion_fraction, "Maximum fraction of deletions for locus to be callable", 0.001, 5));
        addDialogComponent(new DialogComponentNumber(pcr_error, "PCR error rate", 0.0001, 5));
		addDialogComponent(new DialogComponentNumber(call_min_confidence, "Confidence threshold for calling", 1, 5));
        addDialogComponent(new DialogComponentNumber(emit_min_confidence, "Confidence threshold for emitting", 1, 5));
        
        createNewGroup("Malformed read filter");
        addDialogComponent(new DialogComponentBoolean(mbq, "Filter reads with mismatching number of bases and base qualities"));
		        
        //#threads
        createNewGroup("Number of threads");
        addDialogComponent(new DialogComponentNumber(num_threads, "Threads", 1)); 
    }
	
//    private void proxyOptions(){
//	  	createNewTab("Proxy options");
//	  	createNewGroup("General");
//	  	addDialogComponent(new DialogComponentBoolean(useproxy, "Enable proxy"));
//	  	addDialogComponent(new DialogComponentString(proxyhost, "Proxy host"));
//	  	addDialogComponent(new DialogComponentString(proxyport, "Proxy port"));
//	  	createNewGroup("Authentication");
//	  	addDialogComponent(new DialogComponentBoolean(useproxyauth, "Enable authentication"));
//	  	addDialogComponent(new DialogComponentString(proxyuser, "Proxy username"));
//	  	addDialogComponent(new DialogComponentString(proxypassword, "Proxy password"));
//	
//	  	
//	  	useproxy.addChangeListener(new ChangeListener() {
//				public void stateChanged(ChangeEvent e) {
//						if(useproxy.getBooleanValue()){
//							proxyhost.setEnabled(true);
//							proxyport.setEnabled(true);
//							useproxyauth.setEnabled(true);
//						}else{
//							proxyhost.setEnabled(false);
//							proxyport.setEnabled(false);
//							proxyuser.setEnabled(false);
//							proxypassword.setEnabled(false);
//							useproxyauth.setEnabled(false);
//						}
//				}
//			});
//	  	
//	  	useproxyauth.addChangeListener(new ChangeListener() {
//				public void stateChanged(ChangeEvent e) {
//						if(useproxy.getBooleanValue() && useproxyauth.getBooleanValue()){
//							proxyuser.setEnabled(true);
//							proxypassword.setEnabled(true);
//						}else{
//							proxypassword.setEnabled(false);
//							proxyuser.setEnabled(false);
//						}
//				}
//			});
//  }
}

