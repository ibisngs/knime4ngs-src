package de.helmholtz_muenchen.ibis.ngs.gatkbaserecalibration;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "GATKBaseRecalibration" Node.
 * 
 *
 * @author 
 */
public class GATKBaseRecalibrationNodeFactory 
        extends NodeFactory<GATKBaseRecalibrationNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GATKBaseRecalibrationNodeModel createNodeModel() {
        return new GATKBaseRecalibrationNodeModel();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNrNodeViews() {
        return 1;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeView<GATKBaseRecalibrationNodeModel> createNodeView(final int viewIndex,
            final GATKBaseRecalibrationNodeModel nodeModel) {
        return new GATKBaseRecalibrationNodeView(nodeModel);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean hasDialog() {
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeDialogPane createNodeDialogPane() {
        return new GATKBaseRecalibrationNodeDialog();
    }

}

