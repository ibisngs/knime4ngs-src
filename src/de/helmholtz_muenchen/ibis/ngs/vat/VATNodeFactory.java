package de.helmholtz_muenchen.ibis.ngs.vat;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "VAT" Node.
 * 
 *
 * @author 
 */
public class VATNodeFactory 
        extends NodeFactory<VATNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public VATNodeModel createNodeModel() {
        return new VATNodeModel();
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
    public NodeView<VATNodeModel> createNodeView(final int viewIndex,
            final VATNodeModel nodeModel) {
        return new VATNodeView(nodeModel);
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
        return new VATNodeDialog();
    }

}

