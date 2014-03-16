package de.helmholtz_muenchen.ibis.ngs.picard;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentButtonGroup;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/**
 * <code>NodeDialog</code> for the "PicardTools" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class PicardToolsNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring PicardTools node dialog.
     * This is just a suggestion to demonstrate possible default dialog
     * components.
     */
	
	//general options
	final SettingsModelString ptool = new SettingsModelString(PicardToolsNodeModel.CFGKEY_PTOOL, "");	//tool selection
    final SettingsModelString bsformat = new SettingsModelString(PicardToolsNodeModel.CFGKEY_BSFORMAT, PicardToolsNodeModel.DEF_BSFORMAT);	// output format
    final SettingsModelBoolean index = new SettingsModelBoolean(PicardToolsNodeModel.CFGKEY_INDEX, PicardToolsNodeModel.DEF_INDEX);	//create index file?
    final SettingsModelString valstring = new SettingsModelString(PicardToolsNodeModel.CFGKEY_VALSTRING, PicardToolsNodeModel.DEF_VALSTRING);	//validation stringency
	
    //add or replace read groups
    final SettingsModelBoolean use_file_name = new SettingsModelBoolean(PicardToolsNodeModel.CFGKEY_USE_FILE_NAME, PicardToolsNodeModel.DEF_USE_FILE_NAME);	//create id, library, sample name from name of inputfile
    final SettingsModelString id_name = new SettingsModelString(PicardToolsNodeModel.CFGKEY_ID_NAME, PicardToolsNodeModel.DEF_ID_NAME);	//id name
    final SettingsModelString library_name= new SettingsModelString(PicardToolsNodeModel.CFGKEY_LIBRARY_NAME, PicardToolsNodeModel.DEF_LIBRARY_NAME);	//library name
    final SettingsModelString sample_name = new SettingsModelString(PicardToolsNodeModel.CFGKEY_SAMPLE_NAME, PicardToolsNodeModel.DEF_SAMPLE_NAME);	//sample name
    final SettingsModelString platform_unit = new SettingsModelString(PicardToolsNodeModel.CFGKEY_PLATFROM_UNIT, PicardToolsNodeModel.DEF_PLATFORM_UNIT);
    final SettingsModelString platform = new SettingsModelString(PicardToolsNodeModel.CFGKEY_PLATFORM, PicardToolsNodeModel.DEF_PLATFROM);	//sequencing platform
    
    //collect insert size metrics
    final SettingsModelString acc_level = new SettingsModelString(PicardToolsNodeModel.CFGKEY_ACC_LEVEL, PicardToolsNodeModel.DEF_ACC_LEVEL);	//accumulation level
    final SettingsModelBoolean ass_sorted_sm = new SettingsModelBoolean(PicardToolsNodeModel.CFGKEY_ASS_SORTED_SM, PicardToolsNodeModel.DEF_ASS_SORTED_SM);	//input file sorted?
    final SettingsModelDoubleBounded min_pct = new SettingsModelDoubleBounded(PicardToolsNodeModel.CFGKEY_MIN_PCT, PicardToolsNodeModel.DEF_MIN_PCT, PicardToolsNodeModel.MIN_MIN_PC, PicardToolsNodeModel.MAX_MIN_PC); //percentage for discarding reads
    final SettingsModelDoubleBounded deviation = new SettingsModelDoubleBounded(PicardToolsNodeModel.CFGKEY_DEVIATION, PicardToolsNodeModel.DEF_DEVIATION, PicardToolsNodeModel.MIN_DEVIATION, PicardToolsNodeModel.MAX_DEVIATION); //deviation for calculation
    
    //mark duplicates options
    final SettingsModelBoolean remove_dupl = new SettingsModelBoolean(PicardToolsNodeModel.CFGKEY_REMOVE_DUPL, PicardToolsNodeModel.DEF_REMOVE_DUPL);	//remove pcr duplicates
    final SettingsModelBoolean ass_sorted_rd = new SettingsModelBoolean(PicardToolsNodeModel.CFGKEY_ASS_SORTED_RD, PicardToolsNodeModel.DEF_ASS_SORTED_RD);	//input file sorted?
    
    //sorting options
    final SettingsModelString sort_order=new SettingsModelString(PicardToolsNodeModel.CFGKEY_SORT_ORDER, PicardToolsNodeModel.DEF_SORT_ORDER);	//sort order

    
    protected PicardToolsNodeDialog() {
        super();
        
        setUpGeneralOptions();
        setUpReadGroupOptions();
        setUpInsertSizeMetricsOpetions();
        setUpMarkDuplOptions();
        setUpSortingOptions();


        //add Change listener for tool selection -> enables and disables tabs
        ptool.addChangeListener(new ChangeListener(){
        	public void stateChanged(ChangeEvent e){
        		
        		if(ptool.getStringValue().equals("AddOrReplaceReadGroups")){
        			setEnabled(true, "AddOrReplaceReadGroups");
        			bsformat.setEnabled(true);
        			if(bsformat.getStringValue().equals("bam")){
        				index.setEnabled(true);
        			}
    			}
        		else{
        			setEnabled(false, "AddOrReplaceReadGroups");
        		}
        		
        		if(ptool.getStringValue().equals("CollectInsertSizeMetrics")){
        			setEnabled(true, "CollectInsertSizeMetrics");
        			acc_level.setEnabled(true);
        			ass_sorted_sm.setEnabled(true);
        			min_pct.setEnabled(true);
        			deviation.setEnabled(true);
        			bsformat.setEnabled(false);
        			index.setEnabled(false);
        		}
        		else{
        			setEnabled(false, "CollectInsertSizeMetrics");
        			acc_level.setEnabled(false);
        			ass_sorted_sm.setEnabled(false);
        			min_pct.setEnabled(false);
        			deviation.setEnabled(false);
        			bsformat.setEnabled(true);
        			if(bsformat.getStringValue().equals("bam")){
        				index.setEnabled(true);
    				}
        			
        		}
        		
        		if(ptool.getStringValue().equals("MarkDuplicates")){
        			setEnabled(true, "MarkDuplicates");
        			remove_dupl.setEnabled(true);
        			ass_sorted_rd.setEnabled(true);
        			bsformat.setEnabled(true);
        			if(bsformat.getStringValue().equals("bam")){
        				index.setEnabled(true);
        			}
        		}
        		else{
        			setEnabled(false, "MarkDuplicates");
        			remove_dupl.setEnabled(false);
        			ass_sorted_rd.setEnabled(false);
        		}
        		
        		if(ptool.getStringValue().equals("SortSam")){
        			setEnabled(true, "SortSam");
        			sort_order.setEnabled(true);
        			bsformat.setEnabled(true);
        			if(bsformat.getStringValue().equals("bam")){
        				index.setEnabled(true);
        			}
        		}
        		else{
        			setEnabled(false, "SortSam");
        			sort_order.setEnabled(false);
        		}
        		
        	}
        }
        );
    }
    
    private void setUpGeneralOptions(){
    	        
        //tool selection
        createNewGroup("PicardTool selection");
        addDialogComponent(new DialogComponentStringSelection(ptool, "Tool", PicardToolsNodeModel.TOOLS_AVAILABLE));
        //output format
        createNewGroup("Output");
        addDialogComponent(new DialogComponentButtonGroup(bsformat, "Choose output format", true, new String[]{"BAM (required for variant calling)", "SAM"}, new String[]{"bam", "sam"}));
        //create index file?
        addDialogComponent(new DialogComponentBoolean(index, "Create index file (required for GATK)"));
        //validation stringency
        createNewGroup("Validation stringency");
        addDialogComponent(new DialogComponentStringSelection(valstring, "Stringency mode", PicardToolsNodeModel.VALSTRING_AVAILABLE));
        
        // change listener that disables/enables index tick box
        bsformat.addChangeListener(new ChangeListener(){
        	public void stateChanged(ChangeEvent e){
        		if(bsformat.getStringValue().equals("bam")){
        			index.setEnabled(true);
        			index.setBooleanValue(true);
        		}
        		else{
        			index.setEnabled(false);
        			index.setBooleanValue(false);
        		}
        	}
        });
    	
    }
    
    private void setUpMarkDuplOptions(){
        createNewTab("MarkDuplicates");
        setEnabled(false, "MarkDuplicates");
        
        //remove duplicates
        createNewGroup("Removal of PCR duplicates");
        addDialogComponent(new DialogComponentBoolean(remove_dupl, "Do not write duplicated reads to output file"));
        
        //assume sorted
        createNewGroup("Sort order of input file");
        addDialogComponent(new DialogComponentBoolean(ass_sorted_rd, "Input file is sorted by genomic position of the reads"));
    }
    
    private void setUpSortingOptions(){
        createNewTab("SortSam");
        setEnabled(false, "SortSam");
                
        //sort order
        createNewGroup("Sort order");
        addDialogComponent(new DialogComponentButtonGroup(sort_order, "Sort by", true, new String[]{"Mapping position","Read name","Do not sort"}, PicardToolsNodeModel.SORT_AVAILABLE));
        
    }
    
    private void setUpInsertSizeMetricsOpetions(){
        createNewTab("CollectInsertSizeMetrics");
        setEnabled(false, "CollectInsertSizeMetrics");
        
        //metric accumulation level String values
        createNewGroup("Accumulation level for insert size calculation");
        addDialogComponent(new DialogComponentStringSelection(acc_level, "Accumulation level", "ALL_READS", "SAMPLE", "LIBRARY", "READ_GROUP"));
        
        //assume sorted
        createNewGroup("Sort order of input file");
        addDialogComponent(new DialogComponentBoolean(ass_sorted_sm, "Input file is sorted by genomic position of the reads"));

        //minimum pct double zwischen 0 und 1
        createNewGroup("Discard reads of categories RF, TANDEM, FR");
        addDialogComponent(new DialogComponentNumber(min_pct, "Percentage (of overall reads) threshold", 0.01, 3));
        
        //deviation
        createNewGroup("Deviation for plotting and calculation of statistics");
        addDialogComponent(new DialogComponentNumber(deviation, "Deviation", 0.1, 3));
    	
    }

    private void setUpReadGroupOptions(){
        createNewTab("AddOrReplaceReadGroups");
        setEnabled(true, "AddOrReplaceReadGroups");
        
        addDialogComponent(new DialogComponentBoolean(use_file_name, "Use file name to create RG tag"));

        createNewGroup("Read group information");
        
        //ID of RG tag
        addDialogComponent(new DialogComponentString(id_name, "ID name"));
        
        //library name of RG tag
        addDialogComponent(new DialogComponentString(library_name, "Library name"));
        
        //sample name of RG tag
        addDialogComponent(new DialogComponentString(sample_name, "Sample name"));
        
        //sequencing platform unit
        addDialogComponent(new DialogComponentString(platform_unit,"Platform unit"));
        
        //sequencing platform selection
        addDialogComponent(new DialogComponentStringSelection(platform, "Sequencing platform", "ILLUMINA", "SOLID", "LS454", "HELICOS", "PACBIO"));
        
        // change listener to disable/enable name fields
        use_file_name.addChangeListener(new ChangeListener(){
        	public void stateChanged(ChangeEvent e){
        		
        		if(use_file_name.getBooleanValue()){
        			id_name.setEnabled(false);
        			library_name.setEnabled(false);
        			sample_name.setEnabled(false);
        			platform_unit.setEnabled(false);
        		}
        		
        		else{
        			id_name.setEnabled(true);
        			library_name.setEnabled(true);
        			sample_name.setEnabled(true);
        			platform_unit.setEnabled(true);
        		}
        	}
        });
    }

}