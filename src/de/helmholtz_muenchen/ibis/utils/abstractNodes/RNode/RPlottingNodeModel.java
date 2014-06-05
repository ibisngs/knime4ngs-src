package de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

import org.knime.core.data.DataType;
import org.knime.core.data.image.png.PNGImageCell;
import org.knime.core.data.image.png.PNGImageContent;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelColumnName;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.port.PortObject;
import org.knime.core.node.port.PortType;
import org.knime.core.node.port.image.ImagePortObject;
import org.knime.core.node.port.image.ImagePortObjectSpec;

import de.helmholtz_muenchen.ibis.utils.Global;

public abstract class RPlottingNodeModel extends RNodeModel {
	public static final String[] COLOR_SCALES = new String[]{"default", "BrBG (n=11, type=div)", "PiYG (n=11, type=div)", "PRGn (n=11, type=div)", "PuOr (n=11, type=div)", "RdBu (n=11, type=div)", "RdGy (n=11, type=div)", "RdYlBu (n=11, type=div)", "RdYlGn (n=11, type=div)", "Spectral (n=11, type=div)", "Accent (n=8, type=qual)", "Dark2 (n=8, type=qual)", "Paired (n=12, type=qual)", "Pastel1 (n=9, type=qual)", "Pastel2 (n=8, type=qual)", "Set1 (n=9, type=qual)", "Set2 (n=8, type=qual)", "Set3 (n=12, type=qual)", "Blues (n=9, type=seq)", "BuGn (n=9, type=seq)", "BuPu (n=9, type=seq)", "GnBu (n=9, type=seq)", "Greens (n=9, type=seq)", "Greys (n=9, type=seq)", "Oranges (n=9, type=seq)", "OrRd (n=9, type=seq)", "PuBu (n=9, type=seq)", "PuBuGn (n=9, type=seq)", "PuRd (n=9, type=seq)", "Purples (n=9, type=seq)", "RdPu (n=9, type=seq)", "Reds (n=9, type=seq)", "YlGn (n=9, type=seq)", "YlGnBu (n=9, type=seq)", "YlOrBr (n=9, type=seq)", "YlOrRd (n=9, type=seq)"};



	// LOGGER
	private static final NodeLogger LOGGER = NodeLogger.getLogger(RPlottingNodeModel.class);
	
    // cfgkeys
	static public final String CFGKEY_COLUMN_X             = "x column";
	static public final String CFGKEY_COLUMN_Y             = "y column";

	
	static public final String CFGKEY_COLUMN_COLOR         = "color";
	static public final String CFGKEY_LAB_COLOR            = "color label";
	static public final String CFGKEY_SCALE_COLOR          = "color scale";
	static public final String CFGKEY_COLUMN_FILL          = "fill";
	static public final String CFGKEY_LAB_FILL             = "fill label";
	static public final String CFGKEY_SCALE_FILL           = "fill scale";
	static public final String CFGKEY_COLUMN_SHAPE         = "Shape";
	static public final String CFGKEY_LAB_SHAPE            = "Shape label";
	
	static public final String CFGKEY_FACET_GRID           = "facet grid";
	static public final String CFGKEY_FACET_X_COLUMN       = "facet grid (x)";
	static public final String CFGKEY_FACET_Y_COLUMN       = "facet grid (y)";
	
	static public final String CFGKEY_TITLE                = "title";
	static public final String CFGKEY_LAB_X                = "x label";
	static public final String CFGKEY_LAB_Y                = "y label";
	
	static public final String CFGKEY_OUTPUT_PATH          = "output path";
    static public final String CFGKEY_PLOT_WIDTH           = "width";
    static public final String CFGKEY_PLOT_HEIGHT          = "height";
    
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

//	protected RPlottingNodeModel(int nrInDataPorts, int nrOutDataPorts, String script, String[] input_file_arguments, String[] output_file_arguments) {
//		super(nrInDataPorts, nrOutDataPorts, script, input_file_arguments,output_file_arguments);
//	}
//
//	protected RPlottingNodeModel(final PortType[] inPortTypes, final PortType[] outPortTypes, String script, String[] input_file_arguments, String[] output_file_arguments) {
//		super(inPortTypes, outPortTypes, script, input_file_arguments, output_file_arguments);
//	}
	
	protected RPlottingNodeModel(int nrInDataPorts, String script, String[] input_file_arguments, String[] output_file_arguments) {
		super(Global.createOPOs(nrInDataPorts), new PortType[]{ImagePortObject.TYPE}, script, input_file_arguments, output_file_arguments);
	}
	

	@Override
	protected ImagePortObject[] execute(final PortObject[] inObjects, final ExecutionContext exec) throws Exception {
		BufferedDataTable[] inData = new BufferedDataTable[inObjects.length];
		for(int i=0; i<inData.length; i++){
				inData[i] = (BufferedDataTable)inObjects[i];
		}
		
		if(m_col_x.isEnabled()){
			this.addArgument("--col.x", m_col_x);
		}
		if(m_col_y.isEnabled()){
			this.addArgument("--col.y", m_col_y);
		}
		
		if(m_col_color.getColumnName() != null || m_col_color.useRowID()){
			this.addArgument("--col.color", m_col_color);		
			this.addArgument("--lab.color", m_lab_color.getStringValue());
			this.addArgument("--scale.color", m_scale_color.getStringValue().replaceAll( "\\s+.*$", "" ));
		}
		if(m_col_fill.getColumnName() != null || m_col_fill.useRowID()){
			this.addArgument("--col.fill", m_col_fill);		
			this.addArgument("--lab.fill", m_lab_fill.getStringValue());
			this.addArgument("--scale.fill", m_scale_fill.getStringValue().replaceAll( "\\s+.*$", "" ));
		}
		if(m_col_shape.getColumnName() != null || m_col_color.useRowID()){
			this.addArgument("--col.shape", m_col_shape);		
			this.addArgument("--lab.shape", m_lab_shape.getStringValue().replace("\\w.*", ""));
		}
		
		
		if(m_facet.getBooleanValue()){
			this.addArgument("--facet.x", m_facet_xColumn);
			this.addArgument("--facet.y", m_facet_yColumn);
		}
		
		this.addArgument("--title", m_title.getStringValue());
		this.addArgument("--lab.x" , m_lab_x.getStringValue());
		this.addArgument("--lab.y" , m_lab_y.getStringValue());
		
		this.addArgument("--width", m_plotWidth.getIntValue());
		this.addArgument("--height", m_plotHeight.getIntValue());
		this.addArgument("--image", m_output.getStringValue());
		
		
		//////////////////////////////////////////////////////////////////////////
		// PREPARE INPUT DATA
		//////////////////////////////////////////////////////////////////////////
		super.prepareInputData(inData, exec);

		//////////////////////////////////////////////////////////////////////////
		// RUN COMMAND
		//////////////////////////////////////////////////////////////////////////
		exec.setProgress(0.10);
		exec.setProgress("executing script");
		LOGGER.info("Running Rscript with arguments: " + getArgumentsAsVector());
		super.executeScript(exec, null);

		//////////////////////////////////////////////////////////////////////////
		// READ DATA
		//////////////////////////////////////////////////////////////////////////
		exec.setProgress(0.90);
		exec.setProgress("reading output");
		ImagePortObject[] output = null;
		try {
			File imageFile = new File(this.m_output.getStringValue().replace("~",System.getProperty("user.home")));
			FileInputStream imageIS = new FileInputStream(imageFile);
			PNGImageContent imageContent = new PNGImageContent(imageIS);
			ImagePortObjectSpec spec = new ImagePortObjectSpec(DataType.getType( PNGImageCell.class));
			output = new ImagePortObject[]{new ImagePortObject(imageContent, spec)};
		} catch (IOException e) {
			exec.setMessage("Error while reading the image!\n" + e.getMessage());
			throw(e);
		}
		

		exec.setProgress(1.0);

		return(output);

	}


	/////////////////////////////////////////////////////////////////////////////////////////////////
	/// OVERIDE KNIME NODE METHODS
	/////////////////////////////////////////////////////////////////////////////////////////////////
	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void reset() {
		super.reset();
	}

	   /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	super.saveSettingsTo(settings);
    	m_col_x.saveSettingsTo(settings);
    	m_col_y.saveSettingsTo(settings);
    	
    	m_col_color.saveSettingsTo(settings);
    	m_lab_color.saveSettingsTo(settings);
    	m_scale_color.saveSettingsTo(settings);
    	m_col_fill.saveSettingsTo(settings);
    	m_lab_fill.saveSettingsTo(settings);
    	m_scale_fill.saveSettingsTo(settings);
    	m_col_shape.saveSettingsTo(settings);
    	m_lab_shape.saveSettingsTo(settings);
    	
    	m_facet.saveSettingsTo(settings);    	
    	m_facet_xColumn.saveSettingsTo(settings);    	
    	m_facet_yColumn.saveSettingsTo(settings);    	
    	
    	m_title.saveSettingsTo(settings);    	    	    	
    	m_lab_x.saveSettingsTo(settings);
    	m_lab_y.saveSettingsTo(settings);
    	
    	m_plotWidth.saveSettingsTo(settings);
    	m_plotHeight.saveSettingsTo(settings);
    	m_output.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings) throws InvalidSettingsException {
    	super.loadValidatedSettingsFrom(settings);
    	m_col_x.loadSettingsFrom(settings);
    	m_col_y.loadSettingsFrom(settings);
    	
    	m_col_color.loadSettingsFrom(settings);
    	m_lab_color.loadSettingsFrom(settings);
    	m_scale_color.loadSettingsFrom(settings);
    	m_col_fill.loadSettingsFrom(settings);
    	m_lab_fill.loadSettingsFrom(settings);
    	m_scale_fill.loadSettingsFrom(settings);
    	m_col_shape.loadSettingsFrom(settings);
    	m_lab_shape.loadSettingsFrom(settings);
    	
    	m_facet.loadSettingsFrom(settings);
    	m_facet_xColumn.loadSettingsFrom(settings);
    	m_facet_yColumn.loadSettingsFrom(settings);
    	
    	m_title.loadSettingsFrom(settings);
    	m_lab_x.loadSettingsFrom(settings);
    	m_lab_y.loadSettingsFrom(settings);
    	
    	m_plotWidth.loadSettingsFrom(settings);
    	m_plotHeight.loadSettingsFrom(settings);
    	m_output.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings) throws InvalidSettingsException {
    	super.validateSettings(settings);
    	m_col_x.validateSettings(settings);
    	m_col_y.validateSettings(settings);
    	
    	m_col_color.validateSettings(settings);
    	m_lab_color.validateSettings(settings);
    	m_scale_color.validateSettings(settings);
    	m_col_fill.validateSettings(settings);
    	m_lab_fill.validateSettings(settings);
    	m_scale_fill.validateSettings(settings);
    	m_col_shape.validateSettings(settings);
    	m_lab_shape.validateSettings(settings);
    	
    	m_facet.validateSettings(settings);
    	m_facet_xColumn.validateSettings(settings);
    	m_facet_yColumn.validateSettings(settings);
    	
    	m_title.validateSettings(settings);
    	m_lab_x.validateSettings(settings);
    	m_lab_y.validateSettings(settings);
    	
    	m_plotWidth.validateSettings(settings);
    	m_plotHeight.validateSettings(settings);
    	m_output.validateSettings(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
    	super.loadInternals(internDir, exec);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir, final ExecutionMonitor exec) throws IOException, CanceledExecutionException {
    	super.saveInternals(internDir, exec);
    }
	



}
