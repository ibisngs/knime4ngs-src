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
package de.helmholtz_muenchen.ibis.utils.ngs;

import java.util.Arrays;

import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.port.PortType;
import org.knime.core.node.port.PortTypeRegistry;

public class OptionalPorts {

	public static final PortType OPTIONAL_PORT_TYPE = PortTypeRegistry.getInstance().getPortType(BufferedDataTable.class, true);
	
    public static PortType[] createOPOs(final int nrDataPorts, final int... optionalPortsIds)
    {
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
	
}
