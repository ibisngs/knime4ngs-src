package de.helmholtz_muenchen.ibis.utils.ngs;

import java.util.Arrays;

import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.port.PortType;

public class OptionalPorts {

	public static final PortType OPTIONAL_PORT_TYPE = new PortType(BufferedDataTable.class, true);
	
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
