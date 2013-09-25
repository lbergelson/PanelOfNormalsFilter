package org.broadinstitute.cga.tools.gatk.walkers.cancer.filterbypon;

import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFUtils;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: louisb
 * Date: 9/24/13
 * Time: 1:27 PM
 * To change this template use File | Settings | File Templates.
 */
public class FilterByPON extends RodWalker<Integer, Integer>{

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Input(fullName="panel_of_normals", shortName="panel", doc="Panel of normals")
    public List<RodBinding<VariantContext>> panel = Collections.emptyList();;

    @Output(fullName="vcfout", shortName="out", doc="File to write variants which are not filtered out." )
    protected VariantContextWriter vcfWriter = null;




    @Override
    public void initialize(){
        initializeVcfHeader();
    }

    private void initializeVcfHeader() {
        List<String> rodNames = Arrays.asList(variantCollection.variants.getName());
        Map<String, VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(this.getToolkit(), rodNames);

        Set<String> vcfSamples = SampleUtils.getSampleList(vcfRods, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), true);
        vcfWriter.writeHeader(new VCFHeader(headerLines, vcfSamples));
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        //check if we're at the last position
        if ( ref == null){
            return 0;
        }

        final int WIDTH = 5;

        GenomeLoc currentLoc = context.getLocation();
        GenomeLoc window = currentLoc.setStart(currentLoc.setStop(currentLoc,currentLoc.getStop() + WIDTH ), currentLoc.getStart() - WIDTH);

        Collection<VariantContext> panelVariants = tracker.getValues(panel,window);
        Collection<VariantContext> variants = tracker.getValues(variantCollection.variants, currentLoc) ;


        if ( ! reject(variants, panelVariants)) {

            for( VariantContext v : variants ){
                vcfWriter.add(v);
            }

            return 0;
        } else {
            //variant was filtered, so add 1 to the filtered count
            return 1;
        }

    }

    private boolean reject(Collection<VariantContext> variants, Collection<VariantContext> panelVariants) {
        return panelVariants.size() >= 2;
    }

    @Override
    public Integer reduceInit() {
        return 0;
    }

    @Override
    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }
}

