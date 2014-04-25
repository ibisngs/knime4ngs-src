package de.helmholtz_muenchen.ibis.utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataValue;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.port.PortType;




public class Global {
	static public String[] R_OPTIONS_MULTIPLE_TESTING = new String[] {"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"};
	static public String[] R_OPTIONS_GLM_FAMILY       = new String[] {"binomial", "gaussian", "Gamma", "inverse.gaussian", "poisson", "quasi", "quasibinomial", "quasipoisson"};

    static public String FILE_CHOOSER_HISTORIE_ID     = "file.chooser.history.id";
	 public static final PortType OPTIONAL_PORT_TYPE = new PortType(BufferedDataTable.class, true);
	 
	 
	 public static PortType[] createOPOs(final int nrDataPorts, final int... optionalPortsIds){
	     PortType[] portTypes = new PortType[nrDataPorts];
	     Arrays.fill(portTypes, BufferedDataTable.TYPE);        

	     if (optionalPortsIds.length > 0) {
	         for (int portId : optionalPortsIds) {
	             if ((portId - 1) < nrDataPorts) {
	                 portTypes[portId - 1] = OPTIONAL_PORT_TYPE;
	             }
	         }
	     }
	     return portTypes;
	 } 
	 
	 /**
	  * Get the names of columns which are compatible with the specified type
	  * @param specs dataTableSpec
	  * @param type the required Column Type (extends DataValue)
	  * @return
	  */
	 public static String[] getValidCols(DataTableSpec specs, Class<? extends DataValue> type){
		 ArrayList<String> validCols = new ArrayList<String>();
		 Iterator<DataColumnSpec> it = specs.iterator();
		 while(it.hasNext()){
			 DataColumnSpec col = it.next();
			 if(col.getType().isCompatible(type)){
				 validCols.add(col.getName());
			 }
		 }
		 return(validCols.toArray(new String[validCols.size()]));
	 }

	public static ArrayList<String> getVariableTypes(final DataTableSpec inSpec){
		ArrayList<String> classes = new ArrayList<String>(inSpec.getNumColumns());
		for(int c=0; c<inSpec.getNumColumns(); c++){
			classes.add(c, inSpec.getColumnSpec(c).getType().getPreferredValueClass().getName());
		}
		return(classes);
	}

	public static HashSet<String> getUniqueVariableTypes(final DataTableSpec inSpec){
		HashSet<String> classes = new HashSet<String>();
		for(int c=0; c<inSpec.getNumColumns(); c++){
			classes.add(inSpec.getColumnSpec(c).getType().getPreferredValueClass().getName());
		}
		return(classes);
	}
}
