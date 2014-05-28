package de.helmholtz_muenchen.ibis.plotting.histogram;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDouble;

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
public class HistogramNodeDialog extends RPlottingNodeDialog {

	 // models
    protected final SettingsModelBoolean m_density       = new SettingsModelBoolean(HistogramNodeModel.CFGKEY_DENSITY, false);
    protected final SettingsModelBoolean m_density_curve = new SettingsModelBoolean(HistogramNodeModel.CFGKEY_DENSITY_CURVE, false);
    protected final SettingsModelDouble m_binwidth       = new SettingsModelDouble(HistogramNodeModel.CFGKEY_BINWIDTH, 1.0);
   
    /**
     * New pane for configuring the Boxplot node.
     */
	protected HistogramNodeDialog() {
        super();
        // don't need y axis
        d_col_y.getComponentPanel().setVisible(false);
        m_col_y.setEnabled(false);
        
        
        this.createNewTabAt("Histogram Options", 0);
        this.addDialogComponent(new DialogComponentNumber(
        		m_binwidth,
                "binwidth", 1, 5));
        
        this.addDialogComponent(new DialogComponentBoolean(
        		m_density,
                "plot density (instead of counts)"));
        this.addDialogComponent(new DialogComponentBoolean(
        		m_density_curve,
                "add smoothed density line"));
        
        
        m_density.addChangeListener(new ChangeListener() {
      	    public void stateChanged(final ChangeEvent e) {
      	    	m_density_curve.setEnabled(m_density.getBooleanValue());
      	    }
      	});
    }
}

