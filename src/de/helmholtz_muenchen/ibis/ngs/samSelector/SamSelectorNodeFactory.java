package de.helmholtz_muenchen.ibis.ngs.samSelector;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "SamSelector" Node.
 * This Node can be used to select multiple SAM or BAM files.
 *
 * @author Michael Kluge
 */
public class SamSelectorNodeFactory extends NodeFactory<SamSelectorNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public SamSelectorNodeModel createNodeModel() {
        return new SamSelectorNodeModel();
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
    public NodeView<SamSelectorNodeModel> createNodeView(final int viewIndex, final SamSelectorNodeModel nodeModel) {
        return new SamSelectorNodeView(nodeModel);
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
        return new SamSelectorNodeDialog();
    }
}

