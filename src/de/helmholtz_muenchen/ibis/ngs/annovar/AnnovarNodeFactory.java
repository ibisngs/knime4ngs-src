package de.helmholtz_muenchen.ibis.ngs.annovar;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "Annovar" Node.
 * 
 *
 * @author 
 */
public class AnnovarNodeFactory 
        extends NodeFactory<AnnovarNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public AnnovarNodeModel createNodeModel() {
        return new AnnovarNodeModel();
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
    public NodeView<AnnovarNodeModel> createNodeView(final int viewIndex,
            final AnnovarNodeModel nodeModel) {
        return new AnnovarNodeView(nodeModel);
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
        return new AnnovarNodeDialog();
    }

}

