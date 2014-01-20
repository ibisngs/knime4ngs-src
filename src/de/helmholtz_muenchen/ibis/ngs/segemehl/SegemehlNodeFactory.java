package de.helmholtz_muenchen.ibis.ngs.segemehl;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "Segemehl" Node.
 * 
 *
 * @author Jan
 */
public class SegemehlNodeFactory 
        extends NodeFactory<SegemehlNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public SegemehlNodeModel createNodeModel() {
        return new SegemehlNodeModel();
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
    public NodeView<SegemehlNodeModel> createNodeView(final int viewIndex,
            final SegemehlNodeModel nodeModel) {
        return new SegemehlNodeView(nodeModel);
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
        return new SegemehlNodeDialog();
    }

}

