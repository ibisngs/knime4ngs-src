package de.helmholtz_muenchen.ibis.ngs.snpeffgetdb;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "SnpEffGetDB" Node.
 * 
 *
 * @author 
 */
public class SnpEffGetDBNodeFactory 
        extends NodeFactory<SnpEffGetDBNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public SnpEffGetDBNodeModel createNodeModel() {
        return new SnpEffGetDBNodeModel();
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
    public NodeView<SnpEffGetDBNodeModel> createNodeView(final int viewIndex,
            final SnpEffGetDBNodeModel nodeModel) {
        return new SnpEffGetDBNodeView(nodeModel);
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
        return new SnpEffGetDBNodeDialog();
    }

}

