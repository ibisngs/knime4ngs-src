package de.helmholtz_muenchen.ibis.ngs.SNPcall;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "SNPcall" Node.
 * 
 *
 * @author Jan
 */
public class SNPcallNodeFactory 
        extends NodeFactory<SNPcallNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public SNPcallNodeModel createNodeModel() {
        return new SNPcallNodeModel();
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
    public NodeView<SNPcallNodeModel> createNodeView(final int viewIndex,
            final SNPcallNodeModel nodeModel) {
        return new SNPcallNodeView(nodeModel);
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
        return new SNPcallNodeDialog();
    }

}

