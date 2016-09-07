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
import org.knime.core.node.port.PortTypeRegistry;




public class Global {
	static public String[] R_OPTIONS_MULTIPLE_TESTING = new String[] {"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"};
	static public String[] R_OPTIONS_GLM_FAMILY       = new String[] {"binomial", "gaussian", "Gamma", "inverse.gaussian", "poisson", "quasi", "quasibinomial", "quasipoisson"};

    static public String FILE_CHOOSER_HISTORIE_ID     = "file.chooser.history.id";
	 public static final PortType OPTIONAL_PORT_TYPE = PortTypeRegistry.getInstance().getPortType(BufferedDataTable.class, true);
	 
	 
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
