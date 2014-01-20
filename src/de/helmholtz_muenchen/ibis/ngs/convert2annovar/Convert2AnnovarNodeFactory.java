package de.helmholtz_muenchen.ibis.ngs.convert2annovar;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "Convert2Annovar" Node.
 * 
 *
 * @author 
 */
public class Convert2AnnovarNodeFactory 
        extends NodeFactory<Convert2AnnovarNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public Convert2AnnovarNodeModel createNodeModel() {
        return new Convert2AnnovarNodeModel();
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
    public NodeView<Convert2AnnovarNodeModel> createNodeView(final int viewIndex,
            final Convert2AnnovarNodeModel nodeModel) {
        return new Convert2AnnovarNodeView(nodeModel);
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
        return new Convert2AnnovarNodeDialog();
    }

}

