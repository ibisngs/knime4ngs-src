package de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode;

import javax.swing.JFileChooser;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.data.DataValue;
import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentColumnNameSelection;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumberEdit;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelColumnName;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.Global;

/**
 * <code>NodeDialog</code> for the "Boxplot" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Jonas Zierer
 */
public abstract class RPlottingNodeDialog extends DefaultNodeSettingsPane {

    // models
    protected final SettingsModelColumnName m_col_x    	       = new SettingsModelColumnName(RPlottingNodeModel.CFGKEY_COLUMN_X, "");
    protected final SettingsModelColumnName m_col_y    	       = new SettingsModelColumnName(RPlottingNodeModel.CFGKEY_COLUMN_Y, "");
    
    protected final SettingsModelColumnName m_col_color        = new SettingsModelColumnName(RPlottingNodeModel.CFGKEY_COLUMN_COLOR, "");
    protected final SettingsModelString m_lab_color 	       = new SettingsModelString(RPlottingNodeModel.CFGKEY_LAB_COLOR, "");
    protected final SettingsModelString m_scale_color 	       = new SettingsModelString(RPlottingNodeModel.CFGKEY_SCALE_COLOR, RPlottingNodeModel.COLOR_SCALES[0]);
    protected final SettingsModelColumnName m_col_fill         = new SettingsModelColumnName(RPlottingNodeModel.CFGKEY_COLUMN_FILL, "");
    protected final SettingsModelString m_lab_fill  	       = new SettingsModelString(RPlottingNodeModel.CFGKEY_LAB_FILL, "");
    protected final SettingsModelString m_scale_fill 	       = new SettingsModelString(RPlottingNodeModel.CFGKEY_SCALE_FILL, RPlottingNodeModel.COLOR_SCALES[0]);
    protected final SettingsModelColumnName m_col_shape        = new SettingsModelColumnName(RPlottingNodeModel.CFGKEY_COLUMN_SHAPE, "");
    protected final SettingsModelString m_lab_shape 	       = new SettingsModelString(RPlottingNodeModel.CFGKEY_LAB_SHAPE, "");
    
    protected final SettingsModelBoolean m_facet               = new SettingsModelBoolean(RPlottingNodeModel.CFGKEY_FACET_GRID, false);
    protected final SettingsModelColumnName m_facet_xColumn    = new SettingsModelColumnName(RPlottingNodeModel.CFGKEY_FACET_X_COLUMN, "");
    protected final SettingsModelColumnName m_facet_yColumn    = new SettingsModelColumnName(RPlottingNodeModel.CFGKEY_FACET_Y_COLUMN, "");
    
    protected final SettingsModelString m_title   	           = new SettingsModelString(RPlottingNodeModel.CFGKEY_TITLE, "");
    protected final SettingsModelString m_lab_x   	           = new SettingsModelString(RPlottingNodeModel.CFGKEY_LAB_X, "");
    protected final SettingsModelString m_lab_y   	     	   = new SettingsModelString(RPlottingNodeModel.CFGKEY_LAB_Y, "");
	
    protected final SettingsModelInteger m_plotWidth           = new SettingsModelInteger(RPlottingNodeModel.CFGKEY_PLOT_WIDTH, 800);
    protected final SettingsModelInteger m_plotHeight          = new SettingsModelInteger(RPlottingNodeModel.CFGKEY_PLOT_HEIGHT, 600);
    protected final SettingsModelString m_output               = new SettingsModelString(RPlottingNodeModel.CFGKEY_OUTPUT_PATH, "");


    @SuppressWarnings("unchecked")
	protected DialogComponentColumnNameSelection d_col_y = new DialogComponentColumnNameSelection(
  			m_col_y,
  			"Y-axis",
  			0, true, false, DataValue.class);
    /**
     * New pane for configuring the Boxplot node.
     */
	@SuppressWarnings("unchecked")
	protected RPlottingNodeDialog() {
        super();
        
     	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
     	// AXIS SELECTION
     	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        createNewGroup("Column Selection"); 
      	addDialogComponent(new DialogComponentColumnNameSelection(
      			m_col_x,
      			"X-axis",
      			0, true, false, DataValue.class));
      	addDialogComponent(d_col_y);
      	
      	
      	
      	
      	closeCurrentGroup();
        
      	
     	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
     	// LABELS
     	///////////////////////////////////////////////////////////////////////////////////////////////////////////////  	
     	createNewGroup("Labels");
     	addDialogComponent(new DialogComponentString(
     			m_title,
     			"Title Column"));
     	addDialogComponent(new DialogComponentString(
     			m_lab_x,
     			"X Label"));
     	addDialogComponent(new DialogComponentString(
     			m_lab_y,
     			"Y Label"));
     	closeCurrentGroup();
     	
     	
     	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
     	// COLOR
     	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
     	createNewTab("Aesthetics");
     	createNewGroup("Color");
     	m_col_color.addChangeListener(new ChangeListener() {
      	    public void stateChanged(final ChangeEvent e) {
      	    	if(m_col_color.getColumnName() == null && !m_col_color.useRowID()){
      	    		m_lab_color.setEnabled(false);
      	    		m_scale_color.setEnabled(false);
      	    	}else{
      	    		m_lab_color.setEnabled(true);
      	    		m_scale_color.setEnabled(true);
      	    	}
      	    }
      	});
      	addDialogComponent(new DialogComponentColumnNameSelection(
      			m_col_color,
      			"Color Column",
      			0, false, true, DataValue.class));
     	addDialogComponent(new DialogComponentString(
     			m_lab_color,
     			"Color Label"));
     	addDialogComponent(new DialogComponentStringSelection(
     			m_scale_color,
     			"Color Scale",
     			RPlottingNodeModel.COLOR_SCALES));
     	closeCurrentGroup();
     	
     	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
     	// FILL
     	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
     	createNewGroup("Fill Color");
     	m_col_fill.addChangeListener(new ChangeListener() {
      	    public void stateChanged(final ChangeEvent e) {
      	    	if(m_col_fill.getColumnName() == null && !m_col_fill.useRowID()){
      	    		m_lab_fill.setEnabled(false);
      	    		m_scale_fill.setEnabled(false);
      	    	}else{
      	    		m_lab_fill.setEnabled(true);
      	    		m_scale_fill.setEnabled(true);
      	    	}
      	    }
      	});
      	addDialogComponent(new DialogComponentColumnNameSelection(
      			m_col_fill,
      			"Fill Color Column",
      			0, false, true, DataValue.class));
     	addDialogComponent(new DialogComponentString(
     			m_lab_fill,
     			"Fill Color Label"));
     	addDialogComponent(new DialogComponentStringSelection(
     			m_scale_fill,
     			"Fill Color Scale",
     			RPlottingNodeModel.COLOR_SCALES));
     	closeCurrentGroup();
     	
     	
    	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
     	// SHAPE
     	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
     	createNewGroup("Shape");
     	m_col_shape.addChangeListener(new ChangeListener() {
      	    public void stateChanged(final ChangeEvent e) {
      	    	if(m_col_shape.getColumnName() == null && !m_col_shape.useRowID()){
      	    		m_lab_shape.setEnabled(false);
      	    	}else{
      	    		m_lab_shape.setEnabled(true);
      	    	}
      	    }
      	});
      	addDialogComponent(new DialogComponentColumnNameSelection(
      			m_col_shape,
      			"Shape",
      			0, false, true, DataValue.class));
     	addDialogComponent(new DialogComponentString(
     			m_lab_shape,
     			"Shape Label"));
     	closeCurrentGroup();
     	

    	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
     	// OUTPUT
     	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
     	createNewTab("Output Options");
     	addDialogComponent(new DialogComponentFileChooser(
     			m_output, 
     			Global.FILE_CHOOSER_HISTORIE_ID, JFileChooser.SAVE_DIALOG, new String[]{"png"}));
     	addDialogComponent(new DialogComponentNumberEdit(
     			m_plotWidth,
     			"Plot Width"));
     	addDialogComponent(new DialogComponentNumberEdit(
     			m_plotHeight,
     			"Plot Height"));
     	
     	
    	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
     	// FACTES
     	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
     	createNewTab("Facets");
      	m_facet.addChangeListener(new ChangeListener() {
      	    public void stateChanged(final ChangeEvent e) {
      	    	m_facet_xColumn.setEnabled(m_facet.getBooleanValue());
      	    	m_facet_yColumn.setEnabled(m_facet.getBooleanValue());
      	    }
      	});
      	addDialogComponent(new DialogComponentBoolean(
      			m_facet, 
      			"Faceting"));
      	addDialogComponent(new DialogComponentColumnNameSelection(
      			m_facet_xColumn,
      			"X-axis of Facet Grid", 
      			0, false, true, DataValue.class)); 
      	addDialogComponent(new DialogComponentColumnNameSelection(
      			m_facet_yColumn,
      			"Y-axis of Facet Grid", 
      			0, false, true, DataValue.class)); 
   
    }
}

