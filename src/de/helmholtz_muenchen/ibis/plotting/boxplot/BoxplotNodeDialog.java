package de.helmholtz_muenchen.ibis.plotting.boxplot;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.data.DoubleValue;
import org.knime.core.data.IntValue;
import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentButtonGroup;
import org.knime.core.node.defaultnodesettings.DialogComponentColumnFilter;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelFilterString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RPlottingNodeDialog;

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
public class BoxplotNodeDialog extends RPlottingNodeDialog {

    // models
    protected final SettingsModelString m_points                  = new SettingsModelString(BoxplotNodeModel.CFGKEY_POINTS, BoxplotNodeModel.PLOT_POINTS_POSSIBILITIES[0]);
    protected final SettingsModelBoolean m_boxpercol              = new SettingsModelBoolean(BoxplotNodeModel.CFGKEY_BOXPERCOL, false);
    protected final SettingsModelFilterString m_boxpercol_columns = new SettingsModelFilterString(BoxplotNodeModel.CFGKEY_BOXPERCOL_COLUMNS);
   
    /**
     * New pane for configuring the Boxplot node.
     */
	@SuppressWarnings("unchecked")
	protected BoxplotNodeDialog() {
        super();
        this.createNewTabAt("Boxplot Options", 0);
        this.addDialogComponent(new DialogComponentButtonGroup(
        		m_points,
                true, // vertical 
                "plot points",
                BoxplotNodeModel.PLOT_POINTS_POSSIBILITIES));
        m_points.addChangeListener(new ChangeListener() {
      	    public void stateChanged(final ChangeEvent e) {
      	    	if(m_points.getStringValue().equals("all") || m_points.getStringValue().equals("all jittered")){
      	    		m_col_color.setEnabled(true);
      	    		m_lab_color.setEnabled(true);
      	    		m_col_shape.setEnabled(true);
      	    		m_lab_shape.setEnabled(true);
      	    	}else{
      	    		m_col_color.setEnabled(false);
      	    		m_lab_color.setEnabled(false);
      	    		m_col_shape.setEnabled(false);
      	    		m_lab_shape.setEnabled(false);
      	    	}
      	    }
      	});
        
        
        

        this.createNewGroup("Boxes per column");
        addDialogComponent(new DialogComponentBoolean(
        		m_boxpercol, 
        		"Create one box per column"));
        m_boxpercol.addChangeListener(new ChangeListener() {
      	    public void stateChanged(final ChangeEvent e) {
      	    	if(m_boxpercol.getBooleanValue()){
      	    		m_boxpercol_columns.setEnabled(true);
      	    		m_col_x.setEnabled(false);
      	    		m_col_y.setEnabled(false);
      	    	}else{
      	    		m_boxpercol_columns.setEnabled(false);
      	    		m_col_x.setEnabled(true);
      	    		m_col_y.setEnabled(true);
      	    	}
      	    }
      	});
		// Variables
		addDialogComponent(new DialogComponentColumnFilter(
				m_boxpercol_columns,
				0, true, IntValue.class, DoubleValue.class));
		this.closeCurrentGroup();
        
        
    }
}

