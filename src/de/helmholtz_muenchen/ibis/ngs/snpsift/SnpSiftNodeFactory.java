package de.helmholtz_muenchen.ibis.ngs.snpsift;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "SnpSift" Node.
 * 
 *
 * @author Max
 */
public class SnpSiftNodeFactory 
        extends NodeFactory<SnpSiftNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public SnpSiftNodeModel createNodeModel() {
        return new SnpSiftNodeModel();
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
    public NodeView<SnpSiftNodeModel> createNodeView(final int viewIndex,
            final SnpSiftNodeModel nodeModel) {
        return new SnpSiftNodeView(nodeModel);
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
        return new SnpSiftNodeDialog();
    }

}

