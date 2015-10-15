package de.helmholtz_muenchen.ibis.utils.datatypes.file;

import org.knime.core.data.DataType;

@SuppressWarnings("serial")
public class GVCFCell extends FileCell {

	public static final DataType TYPE = DataType.getType(GVCFCell.class);
	
	GVCFCell(String str) {
		super(str);
	}
}
