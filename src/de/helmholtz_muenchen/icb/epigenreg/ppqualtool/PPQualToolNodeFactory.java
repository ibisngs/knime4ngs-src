package de.helmholtz_muenchen.icb.epigenreg.ppqualtool;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "PPQualTool" Node.
 * 
 *
 * @author 
 */
public class PPQualToolNodeFactory 
        extends NodeFactory<PPQualToolNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public PPQualToolNodeModel createNodeModel() {
        return new PPQualToolNodeModel();
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
    public NodeView<PPQualToolNodeModel> createNodeView(final int viewIndex,
            final PPQualToolNodeModel nodeModel) {
        return new PPQualToolNodeView(nodeModel);
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
        return new PPQualToolNodeDialog();
    }

}

