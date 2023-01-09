from Bio import Entrez
import myvariant
class VarintInfoIO:
    def parseMyVarinats(self, request,fabricatedTerm):
        mv = myvariant.MyVariantInfo()
        mvVar = ''
        Entrez.email = "usman.athar@gmail.com"
        handle = Entrez.esearch(db="snp", term=fabricatedTerm, retmax=1)
        variantData = Entrez.read(handle)
        totalSNPs = variantData['Count']
        varids = variantData["IdList"]
        chekboxValues = str(request.POST.getlist('chekboxValues', ''))
        rsIds = []
        var = ''
        print('chekboxValues=== ', chekboxValues)
        fabricatedFields=[]
        for varid in varids:
            id = 'rs' + str(varid)
            var = mv.getvariant(id)
            if chekboxValues.__contains__('All'):
                print('fld=== ', mv.getvariant(id, fields=['dbnsfp.cadd','dbnsfp.cadd.pred', 'dbnsfp.cadd.raw_rankscore',
                                                           'dbnsfp.cadd.raw_score',
                                                           'dbnsfp.sift','dbnsfp.sift.pred',
                                                           'dbnsfp.sift.converted_rankscore', 'dbnsfp.sift.score',
                                                           'dbnsfp.sift4g,dbnsfp.sift4g.pred',
                                                           'dbnsfp.sift4g.converted_rankscore', 'dbnsfp.sift4g.score',
                                                           'dbnsfp.provean,dbnsfp.provean.pred',
                                                           'dbnsfp.provean.rankscore', 'dbnsfp.provean.score',
                                                           'dbnsfp.primateai,dbnsfp.primateai.pred',
                                                           'dbnsfp.primateai.rankscore', 'dbnsfp.primateai.score',
                                                           'dbnsfp.polyphen2.hdiv,dbnsfp.polyphen2.hdiv.pred',
                                                           'dbnsfp.polyphen2.hdiv.rankscore',
                                                           'dbnsfp.polyphen2.hdiv.score',
                                                           'dbnsfp.polyphen2.hvar', 'dbnsfp.polyphen2.hvar.pred',
                                                           'dbnsfp.polyphen2.hvar.rankscore',
                                                           'dbnsfp.polyphen2.hvar.score',
                                                           'dbnsfp.revel', 'dbnsfp.revel.rankscore',
                                                           'dbnsfp.revel.score',
                                                           'dbnsfp.mvp', 'dbnsfp.mvp.rankscore', 'dbnsfp.mvp.score',
                                                           'dbnsfp.mpc', 'dbnsfp.mpc.rankscore', 'dbnsfp.mpc.score',
                                                           'dbnsfp.mutpred,dbnsfp.mutpred.pred',
                                                           'dbnsfp.mutpred.rankscore', 'dbnsfp.mutpred.score',
                                                           'dbnsfp.mutationtaster,dbnsfp.mutationtaster.pred',
                                                           'dbnsfp.mutationtaster.converted_rankscore',
                                                           'dbnsfp.mutationtaster.score',
                                                           'dbnsfp.mutationassessor,dbnsfp.mutationassessor.pred',
                                                           'dbnsfp.mutationassessor.rankscore',
                                                           'dbnsfp.mutationassessor.score',
                                                           'dbnsfp.metarnn,dbnsfp.metarnn.pred',
                                                           'dbnsfp.metarnn.rankscore', 'dbnsfp.metarnn.score',
                                                           'dbnsfp.metasvm,dbnsfp.metasvm.pred',
                                                           'dbnsfp.metasvm.rankscore', 'dbnsfp.metasvm.score',
                                                           'dbnsfp.metalr,dbnsfp.metalr.pred',
                                                           'dbnsfp.metalr.rankscore', 'dbnsfp.metalr.score',
                                                           'dbnsfp.m-cap,dbnsfp.m-cap.pred', 'dbnsfp.m-cap.rankscore',
                                                           'dbnsfp.m-cap.score',
                                                           'dbnsfp.lrt,dbnsfp.lrt.pred',
                                                           'dbnsfp.lrt.converted_rankscore', 'dbnsfp.lrt.score',
                                                           'dbnsfp.list-s2,dbnsfp.list-s2.pred',
                                                           'dbnsfp.list-s2.rankscore', 'dbnsfp.list-s2.score',

                                                           'dbnsfp.fathmm,dbnsfp.fathmm.pred',
                                                           'dbnsfp.fathmm.rankscore', 'dbnsfp.fathmm.score',
                                                           'dbnsfp.fathmm-mkl,dbnsfp.fathmm-mkl.coding_pred',
                                                           'dbnsfp.fathmm-mkl.coding_rankscore',
                                                           'dbnsfp.fathmm-mkl.coding_score',
                                                           'dbnsfp.fathmm-xf,dbnsfp.fathmm-xf.coding_pred',
                                                           'dbnsfp.fathmm-xf.coding_rankscore',
                                                           'dbnsfp.fathmm-xf.coding_score',
                                                           'dbnsfp.bayesdel.add_af,dbnsfp.bayesdel.add_af.pred',
                                                           'dbnsfp.bayesdel.add_af.rankscore',
                                                           'dbnsfp.bayesdel.add_af.score',
                                                           'dbnsfp.bayesdel.no_af,dbnsfp.bayesdel.no_af.pred',
                                                           'dbnsfp.bayesdel.no_af.rankscore',
                                                           'dbnsfp.bayesdel.no_af.score',
                                                           'dbnsfp.aloft,dbnsfp.aloft.pred', 'dbnsfp.aloft.confidence',
                                                           'dbnsfp.vest4', 'dbnsfp.vest4.rankscore',
                                                           'dbnsfp.vest4.score',
                                                           'dbnsfp.dann', 'dbnsfp.dann.rankscore', 'dbnsfp.dann.score',
                                                           'dbnsfp.eigen', 'dbnsfp.eigen.raw_coding',
                                                           'dbnsfp.eigen.raw_coding_rankscore',
                                                           'dbnsfp.eigen-pc', 'dbnsfp.eigen-pc.raw_coding',
                                                           'dbnsfp.eigen-pc.raw_coding_rankscore',
                                                           'dbnsfp.deogen2', 'dbnsfp.deogen2.pred',
                                                           'dbnsfp.deogen2.rankscore', 'dbnsfp.deogen2.score',
                                                           'dbnsfp.genocanyon', 'dbnsfp.genocanyon.rankscore',
                                                           'dbnsfp.genocanyon.score', ]))

            if chekboxValues.__contains__('sift'):
                fabricatedFields.append('dbnsfp.sift.pred, dbnsfp.sift.converted_rankscore, dbnsfp.sift.score')
            if chekboxValues.__contains__('sift4g'):
                fabricatedFields.append('dbnsfp.sift4g.pred, dbnsfp.sift4g.converted_rankscore, dbnsfp.sift4g.score')
            if chekboxValues.__contains__('provean'):
                fabricatedFields.append('dbnsfp.provean.pred,dbnsfp.provean.rankscore, dbnsfp.provean.score')
            if chekboxValues.__contains__('cadd'):
                fabricatedFields.append('dbnsfp.cadd.pred, dbnsfp.cadd.raw_rankscore, dbnsfp.cadd.raw_score')
            if chekboxValues.__contains__('mvp'):
                fabricatedFields.append('dbnsfp.mvp.rankscore', 'dbnsfp.mvp.score')
            if chekboxValues.__contains__('polyphen2hvar'):
                fabricatedFields.append('dbnsfp.polyphen2.hvar.pred, dbnsfp.polyphen2.hvar.rankscore','dbnsfp.polyphen2.hvar.score')
            if chekboxValues.__contains__('polyphen2hdiv'):
                fabricatedFields.append('dbnsfp.polyphen2.hdiv.pred,dbnsfp.polyphen2.hdiv.rankscore,dbnsfp.polyphen2.hdiv.score')
            if chekboxValues.__contains__('primateai'):
                fabricatedFields.append('dbnsfp.primateai.pred, dbnsfp.primateai.rankscore, dbnsfp.primateai.score')
            if chekboxValues.__contains__('revel'):
                fabricatedFields.append('dbnsfp.revel.rankscore, dbnsfp.revel.score')
            if chekboxValues.__contains__('mpc'):
                fabricatedFields.append('dbnsfp.mpc.rankscore, dbnsfp.mpc.score')
            if chekboxValues.__contains__('mutpred'):
                fabricatedFields.append('dbnsfp.mutpred.pred, dbnsfp.mutpred.rankscore, dbnsfp.mutpred.score')
            if chekboxValues.__contains__('mtaster'):
                fabricatedFields.append('dbnsfp.mutationtaster.pred,dbnsfp.mutationtaster.converted_rankscore,dbnsfp.mutationtaster.score')
            if chekboxValues.__contains__('massessor'):
                fabricatedFields.append('dbnsfp.mutationassessor.pred,dbnsfp.mutationassessor.rankscore,dbnsfp.mutationassessor.score')
            if chekboxValues.__contains__('mrnn'):
                fabricatedFields.append('dbnsfp.metarnn.pred, dbnsfp.metarnn.rankscore, dbnsfp.metarnn.score')
            if chekboxValues.__contains__('msvm'):
                fabricatedFields.append('dbnsfp.metasvm.pred, dbnsfp.metasvm.rankscore, dbnsfp.metasvm.score')
            if chekboxValues.__contains__('mlr'):
                fabricatedFields.append('dbnsfp.metalr.pred, dbnsfp.metalr.rankscore, dbnsfp.metalr.score')
            if chekboxValues.__contains__('mcap'):
                fabricatedFields.append('dbnsfp.m-cap.pred, dbnsfp.m-cap.rankscore, dbnsfp.m-cap.score')
            if chekboxValues.__contains__('lrt'):
                fabricatedFields.append('dbnsfp.lrt.pred,dbnsfp.lrt.converted_rankscore, dbnsfp.lrt.score')
            if chekboxValues.__contains__('ls2'):
                fabricatedFields.append('dbnsfp.list-s2.pred, dbnsfp.list-s2.rankscore, dbnsfp.list-s2.score')
            if chekboxValues.__contains__('fathmm'):
                fabricatedFields.append('dbnsfp.fathmm.pred, dbnsfp.fathmm.rankscore, dbnsfp.fathmm.score')
            if chekboxValues.__contains__('fxf'):
                fabricatedFields.append('dbnsfp.fathmm-xf.coding_pred,dbnsfp.fathmm-xf.coding_rankscore, dbnsfp.fathmm-xf.coding_score')
            if chekboxValues.__contains__('fmkl'):
                fabricatedFields.append('fathmm-mkl.coding_pred,dbnsfp.fathmm-mkl.coding_rankscore,dbnsfp.fathmm-mkl.coding_score')
            if chekboxValues.__contains__('bdeladdaf'):
                fabricatedFields.append('dbnsfp.bayesdel.add_af.pred, dbnsfp.bayesdel.add_af.rankscore,dbnsfp.bayesdel.add_af.score')
            if chekboxValues.__contains__('bdelnoaf'):
                fabricatedFields.append('dbnsfp.bayesdel.no_af.pred, dbnsfp.bayesdel.no_af.rankscore,dbnsfp.bayesdel.no_af.score')
            if chekboxValues.__contains__('aloft'):
                fabricatedFields.append('dbnsfp.aloft.pred', 'dbnsfp.aloft.confidence')
            if chekboxValues.__contains__('vest4'):
                fabricatedFields.append('dbnsfp.vest4.rankscore,dbnsfp.vest4.score')
            if chekboxValues.__contains__('dann'):
                fabricatedFields.append('dbnsfp.dann.rankscore', 'dbnsfp.dann.score')
            if chekboxValues.__contains__('eigen'):
                fabricatedFields.append('dbnsfp.eigen.raw_coding, dbnsfp.eigen.raw_coding_rankscore')
            if chekboxValues.__contains__('eigenpc'):
                fabricatedFields.append('dbnsfp.eigen-pc.raw_coding, dbnsfp.eigen-pc.raw_coding_rankscore')
            if chekboxValues.__contains__('doegen2'):
                fabricatedFields.append('dbnsfp.deogen2.pred,dbnsfp.deogen2.rankscore, dbnsfp.deogen2.score')
            if chekboxValues.__contains__('genocanyon'):
                fabricatedFields.append('dbnsfp.genocanyon.rankscore,dbnsfp.genocanyon.score')




            if chekboxValues.__contains__('polyphen'):  # writing GRCh38 chromosome coordinates in csv format
                mv.getvariant(id, fields=['cadd.sift.cat', 'cadd.sift.val'])
            if chekboxValues.__contains__('provean'):
                mv.getvariant(id, fields=['cadd.sift.cat', 'cadd.sift.val'])

                id = 'rs' + str(varid)

                mv = myvariant.MyVariantInfo()
                mvVar = mv.getvariant(id, fields=['cadd.sift.cat', 'cadd.sift.val'])
                # siftFields=mvVar.get_fields("sift")

            print('sift--- ', mv.getvariant(id, fields=fabricatedFields))



