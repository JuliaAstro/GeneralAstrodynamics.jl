"""
Unit tests for various solvers.
"""
module AstrodynamicalSolversTests

using Test, AstrodynamicalSolvers
using AstrodynamicalCalculations
using StaticArrays
using OrdinaryDiffEqVerner
using AstrodynamicalModels
using LinearAlgebra

@testset "Orbit Propagation" begin

    orbit = Orbit(rand(CartesianState), rand(R2BParameters))
    @test propagate(orbit, 1.0) isa ODESolution

    state = copy(orbit.state)
    @test isnothing(propagate!(orbit, 100.0))
    @test orbit.state != state

end

@testset "Lyapunov Orbit Correction" begin

    μ = 0.012150584395829193
    u = richardson_ic(μ, 1)
    u = AstrodynamicalSolvers.CR3BSolvers.lyapunov(u.x, u.ẏ, μ, u.Δt)


    @test u.x ≈ 0.8222791798525636
    @test u.ẏ ≈ 0.13799313228400178
    @test u.Δt ≈ 2.7536820160579087

end

@testset "Halo Orbit Correction" begin

    μ = 0.012150584395829193
    u = richardson_ic(μ, 2; Z = 0.005)
    u = AstrodynamicalSolvers.CR3BSolvers.halo(u.x, u.z, u.ẏ, μ, u.Δt)

    @test u.x ≈ 1.1202340567932783
    @test u.z ≈ 0.004589679675825104
    @test u.ẏ ≈ 0.17648270824601714
    @test u.Δt ≈ 3.4152029027146815

end

@testset "Semantic Orbit Correction" begin

    μ = 0.012150584395829193
    u = halo(μ, 2)

    @test u.x ≈ 1.124357139749168
    @test u.ẏ ≈ 0.15714566174026515
    @test u.Δt ≈ 3.4068306866985814

    u = halo(μ, 1; amplitude = 0.005)
    @test u.x ≈ 0.823388563881332
    @test u.z ≈ 0.005553604696093592
    @test u.ẏ ≈ 0.12683910108732768
    @test u.Δt ≈ 2.7432058155507653

end

@testset "Monodromy Matrices" begin

    μ = 0.012150584395829193
    ic = halo(μ, 1; amplitude = 0.005)
    u = CartesianState(ic)

    @test monodromy(u, μ, ic.Δt, CR3BFunction(stm = true)) ≈ [
        1317.472125300217 -339.26725145920585 -22.23471304682866 388.2455372345198 126.49047353668445 -3.3779594177227
        -429.7764385389933 111.53280315746214 7.269431052432993 -126.49047353678779 -41.404653982215095 1.0977750840404263
        -11.440384647744368 2.9348175664641527 1.1929568568249131 -3.3779594177162653 -1.0977750840374354 -0.05845992955397675
        3612.12986893901 -929.3271447079647 -60.998106970610365 1064.4911782262984 346.9671305741014 -9.244834479682378
        -1482.5514995781887 382.0700769975138 25.039782109090023 -437.2238230101127 -141.4481439160884 3.8211012689780377
        -75.5369690753493 19.429643984545518 1.2797982252155324 -22.23471304678537 -7.269431052412963 1.1929568568243507
    ]
end

@testset "Manifold Computations" begin

    μ = 3.0034805945423635e-6
    u = halo(μ, 1; amplitude = 0.005)

    @test u.x ≈ 0.9894058673149723
    @test u.z ≈ 0.005986079972943528
    @test u.ẏ ≈ 0.01250973206745133
    @test u.Δt ≈ 3.0247281222825437

    T = u.Δt
    u = CartesianState(u)
    Φ = monodromy(u, μ, T, CR3BFunction(stm = true))

    @test Φ ≈ [
        363.5695646260086 -150.84678203289957 -110.89785303818041 147.89291423249102 39.776655263625734 -19.506758521645583
        -98.06514156622599 41.30059843066769 30.40698418319986 -39.776655263719505 -11.022980358853847 5.037110248654976
        -47.32080144069423 19.504183945583154 15.347602546829501 -19.506758521547493 -5.037110248617916 2.2693997822749954
        697.403493226936 -288.5621336671776 -212.7440873239205 284.01625409850567 76.01918084826457 -37.246580943614184
        -355.9757987239261 148.34069711110283 109.04417094750771 -144.93904643230533 -38.25271209672777 19.509333097636393
        -273.5580556895539 112.7515351286987 83.74272666784147 -110.89785303784609 -30.406984183046543 15.347602546859227
    ]

    ics = let
        problem = ODEProblem(CR3BFunction(stm = true), vcat(u, vec(I(6))), (0, T), (μ,))
        solution =
            solve(problem, Vern9(), reltol = 1e-12, abstol = 1e-12, saveat = (T / 10))

        solution.u
    end

    @test ics ≈ [
        [
            0.9894058673149723,
            0.0,
            0.005986079972943528,
            0.0,
            0.01250973206745133,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
        ],
        [
            0.9898744278702267,
            0.0035899577966164516,
            0.005258145836456971,
            0.0029834070095247816,
            0.010601176566460811,
            -0.004737084798616238,
            1.2371458411768335,
            -0.0695323761780384,
            -0.09847482817935505,
            1.5670044822929252,
            -0.6978155777157768,
            -0.6591788390257907,
            -0.035434947222201855,
            0.9316804373498873,
            0.011356684405103354,
            -0.3467997416963855,
            -0.3896842189464582,
            0.11054943644691602,
            -0.09395345250960796,
            0.030697477319750244,
            0.9321583672024584,
            -0.6007887072141127,
            0.3009033764290183,
            -0.4464993538639933,
            0.308242649904382,
            -0.09438488946911146,
            -0.009837529520694738,
            1.0591107163922577,
            -0.6442600720711773,
            -0.09773198568283194,
            0.088131316473712,
            0.2767705969926594,
            0.00023169992649377062,
            0.5611428035212996,
            0.7487829717456858,
            0.0027324590487216266,
            -0.009332241416375939,
            0.0031441665203776454,
            0.295364864880218,
            -0.08956991747299407,
            0.040913095200603176,
            0.9279360834253226,
        ],
        [
            0.9910796788130449,
            0.006066347272099451,
            0.0032124402157081276,
            0.0046840867671740635,
            0.005271523308608312,
            -0.008559051744536353,
            1.950040239465447,
            -0.5899540353747552,
            -0.40535975374341493,
            3.164386111879475,
            -3.0927704802812737,
            -1.3725873609363672,
            -0.26843184364042577,
            0.8374953199128866,
            0.08343688204016234,
            -1.2726857507607194,
            -0.053624055517039464,
            0.3808968730554891,
            -0.34391493779148813,
            0.2355718050038613,
            0.7337936721629754,
            -1.018572831318294,
            1.1365901126247782,
            -0.8528469298727669,
            0.6581421492414853,
            -0.4199981787375973,
            -0.07819233335820817,
            1.296363278219651,
            -1.6024827121626852,
            -0.3800784071360596,
            0.3163004973499631,
            0.4108044051437538,
            0.0016055058382268882,
            0.898197886993142,
            0.07948286238182267,
            0.002678040639188702,
            -0.06481411101658309,
            0.04588962838569348,
            0.5426402313760416,
            -0.2821160224582451,
            0.28560464167961974,
            0.667837207225828,
        ],
        [
            0.9924953609558557,
            0.006527893673766641,
            0.00027703311101019744,
            0.004308405329300428,
            -0.0025418931279188855,
            -0.010403893071012492,
            3.194979782068513,
            -2.2480244592658707,
            -0.8912901561737726,
            5.260531462882106,
            -8.629748694340314,
            -1.6131421475527377,
            -0.8480334714756765,
            1.051120901341595,
            0.23328760581650157,
            -2.6444533370499497,
            1.8137063921561256,
            0.5366273688482629,
            -0.7042926255865667,
            0.7799660829349495,
            0.4234742102702738,
            -1.4067052181979423,
            2.6169072350232443,
            -1.2068964645667555,
            1.1360123079440265,
            -1.169120882563875,
            -0.23731765154027032,
            1.9926119503214768,
            -3.6422821167111397,
            -0.5954430923230926,
            0.6017020436501692,
            0.2894337604412833,
            -0.0024364332298012063,
            0.9513391643697624,
            -0.9600543630574989,
            -0.031952547093191674,
            -0.17594483247914564,
            0.20206265728008083,
            0.6681165100769222,
            -0.4390532271528181,
            0.78324225648997,
            0.09362900661603259,
        ],
        [
            0.9934985373534506,
            0.004414493412947694,
            -0.002691682347757935,
            0.0021433455930179954,
            -0.011393931006874056,
            -0.008319425757109752,
            5.490003988184001,
            -6.487014689628338,
            -0.8875841662056261,
            11.302082651153246,
            -20.6917674268284,
            3.631168290269304,
            -1.9880804273779376,
            2.231086179719102,
            0.24089957990992225,
            -5.357497029086754,
            6.581690137390508,
            -1.1855913442432098,
            -1.309840956076913,
            1.977719338435574,
            -0.07622512726413008,
            -3.064737419133734,
            5.656005016906049,
            -2.517216773282615,
            2.0428301884897455,
            -2.8931589692820436,
            -0.24254153002609002,
            4.620462672370931,
            -8.271485402807551,
            1.3685197351272902,
            0.9017837097616304,
            -0.25264672130310206,
            0.015551302170489048,
            1.1676475337624053,
            -2.821743441950994,
            0.35282832145644405,
            -0.32809129465828046,
            0.5293704124595559,
            0.5416243419739394,
            -0.5956212745474034,
            1.358220989856255,
            -1.0608489042619311,
        ],
        [
            0.9938018000948632,
            3.3814004664628033e-16,
            -0.004091056650145035,
            -1.5713283266960783e-15,
            -0.01638585277717688,
            -1.3586014618110917e-15,
            11.747249011272089,
            -14.478469379121254,
            3.5877588628656816,
            33.92593778407872,
            -27.57659845017894,
            29.9164518117297,
            -4.718395501501515,
            5.047894587113171,
            -1.3552904353555262,
            -14.253386838881683,
            10.544629652277052,
            -11.079226610892352,
            -3.1371461505739706,
            4.163339859874222,
            -1.6641718993686987,
            -10.324059030179551,
            7.616817895784776,
            -9.057160270532727,
            4.635450259594342,
            -6.044958610263879,
            1.5056456364454869,
            14.029780852500256,
            -10.722766001172136,
            11.748164468177174,
            1.536462255983425,
            -1.421038054235922,
            0.5070088326456889,
            3.600270813484911,
            -4.334194398507778,
            3.495249696857971,
            -0.6039468437963879,
            0.9459379702817146,
            -0.05896360053465793,
            -1.4370597187069243,
            1.110717948772804,
            -2.9657277873069456,
        ],
        [
            0.9934985373534488,
            -0.004414493412947038,
            -0.002691682347758814,
            -0.002143345593023276,
            -0.011393931006874526,
            0.008319425757106409,
            27.57759248770883,
            -19.172242333275435,
            16.311421629570358,
            72.18267416063888,
            -0.5217071771052122,
            47.17874452640869,
            -11.340443337730555,
            7.047344305454021,
            -6.297176883156646,
            -30.26726229627077,
            1.2161844662072916,
            -18.973437822540653,
            -8.071483120027205,
            5.474388542930831,
            -5.324321068514165,
            -22.672475850776927,
            0.20393609974423393,
            -13.184980940723491,
            11.118981342445453,
            -7.811335199137715,
            6.498398319499484,
            29.399620552354797,
            0.14240562800692375,
            18.498373101895407,
            3.2910208635699147,
            -2.296036749641137,
            2.046348060177439,
            8.079637992426136,
            -0.9590910848556188,
            5.741127700197823,
            -1.3335010300736867,
            1.0231776131283727,
            -1.1156153587866342,
            -3.5820118833628345,
            -0.6679709374220457,
            -3.564673805710898,
        ],
        [
            0.9924953609558516,
            -0.006527893673766245,
            0.00027703311100841577,
            -0.004308405329310747,
            -0.0025418931279196384,
            0.01040389307101022,
            57.181133640704,
            -16.449867036121947,
            27.93439296013566,
            128.6079397443349,
            12.671869857737473,
            24.502246328517415,
            -23.759695823973058,
            6.304150081118974,
            -11.100619162382069,
            -53.86522120299188,
            -3.8482982939686394,
            -10.515733266818891,
            -17.31335910736424,
            4.773489630454567,
            -8.455598228679394,
            -39.86056875482027,
            -3.118645953595972,
            -6.051829092328109,
            23.153417082539452,
            -6.635970433902056,
            11.051184168659917,
            52.26235486966227,
            5.293684620687887,
            9.554469389008181,
            6.559940802699193,
            -2.1869174797907305,
            3.4351989345430938,
            14.03599732681375,
            0.9722818404698685,
            2.8365521687132196,
            -2.8986711961763487,
            0.6731600249549142,
            -1.968672216943043,
            -7.028498918954745,
            -1.2720733545245055,
            -1.7866195396268072,
        ],
        [
            0.9910796788130367,
            -0.006066347272098975,
            0.003212440215706125,
            -0.004684086767192054,
            0.005271523308610293,
            0.008559051744537721,
            109.60942630891857,
            -15.717566986242204,
            28.48576764394318,
            227.13204414846317,
            -16.067706286794227,
            -25.486055887477125,
            -45.63564975599253,
            6.432439099160851,
            -11.495435787072722,
            -94.48845805804095,
            8.089793565062266,
            9.893609778856955,
            -33.43705608933802,
            4.824834637056938,
            -8.223622040706017,
            -69.50428987878391,
            5.887194656137066,
            8.94800729228516,
            44.47337250008668,
            -6.301270722022312,
            11.207898598128924,
            92.4220929447515,
            -6.419105350976998,
            -10.413143142373471,
            12.237916296308878,
            -2.1875950941690143,
            3.4824234136026435,
            24.51470799457985,
            -1.8729279460958983,
            -2.995290734346036,
            -5.788225941827381,
            0.5091307093232996,
            -2.0570560656612114,
            -12.510339847465445,
            0.6597236178178729,
            1.4334143722627077,
        ],
        [
            0.9898744278702115,
            -0.003589957796614368,
            0.005258145836456293,
            -0.0029834070095553396,
            0.010601176566470631,
            0.004737084798624322,
            201.86456801433138,
            -32.79371839578189,
            8.778449492317803,
            398.7473687287775,
            -112.20521765826136,
            -113.34271959731471,
            -83.91562772874882,
            13.925566566146511,
            -3.5406485779734918,
            -165.18626828653407,
            47.779213285258685,
            46.26199010349747,
            -61.59279313916355,
            10.310376359525133,
            -1.9205492291693678,
            -121.55675882377591,
            35.015239886510024,
            35.31449273448445,
            82.03383012577265,
            -13.217635207287676,
            3.220378725418548,
            162.39550677996976,
            -45.58234688180706,
            -45.886714273972885,
            22.204797865448924,
            -4.032943867089017,
            1.2329333115454304,
            43.16784495801038,
            -11.990747315552705,
            -12.788440637646016,
            -10.823837837259962,
            1.4311210347644943,
            -0.9158767631079093,
            -21.58671912250025,
            6.241725963443699,
            6.540160859407583,
        ],
        [
            0.9894058673149445,
            7.496100275905479e-15,
            0.005986079972946951,
            -5.224004126394785e-14,
            0.012509732067479662,
            2.0208745832464827e-14,
            363.5695646260086,
            -98.06514156622599,
            -47.32080144069423,
            697.403493226936,
            -355.9757987239261,
            -273.5580556895539,
            -150.84678203289957,
            41.30059843066769,
            19.504183945583154,
            -288.5621336671776,
            148.34069711110283,
            112.7515351286987,
            -110.89785303818041,
            30.40698418319986,
            15.347602546829501,
            -212.7440873239205,
            109.04417094750771,
            83.74272666784147,
            147.89291423249102,
            -39.776655263719505,
            -19.506758521547493,
            284.01625409850567,
            -144.93904643230533,
            -110.89785303784609,
            39.776655263625734,
            -11.022980358853847,
            -5.037110248617916,
            76.01918084826457,
            -38.25271209672777,
            -30.406984183046543,
            -19.506758521645583,
            5.037110248654976,
            2.2693997822749954,
            -37.246580943614184,
            19.509333097636393,
            15.347602546859227,
        ],
    ]

    unstable = [diverge(ic[1:6], reshape(ic[7:end], 6, 6), Φ; eps = 1e-8) for ic in ics]

    @test unstable ≈ [
        [
            0.9894058633306991,
            1.072047306563596e-9,
            0.0059860804945415895,
            -7.645981012482653e-9,
            0.012509735977210676,
            2.9928203953031123e-9,
        ],
        [
            0.989874424256154,
            0.0035899593869439046,
            0.0052581467844121385,
            0.0029834002996964354,
            0.010601181942007543,
            -0.004737081701453126,
        ],
        [
            0.9910796756645479,
            0.00606634941623854,
            0.003212441327024458,
            0.004684081166578108,
            0.005271530075233423,
            -0.00855904907981877,
        ],
        [
            0.9924953583827216,
            0.0065278962838587464,
            0.0002770340761658626,
            0.004308400722444655,
            -0.0025418852194094683,
            -0.010403891704544876,
        ],
        [
            0.9934985353165007,
            0.004414496217085891,
            -0.002691681954480171,
            0.002143341012609228,
            -0.011393922984824075,
            -0.008319427338105041,
        ],
        [
            0.9938017980401521,
            2.6371735230836064e-9,
            -0.004091057288364902,
            -6.142464108479443e-9,
            -0.01638584808508263,
            -5.3544911573575455e-9,
        ],
        [
            0.993498534435403,
            -0.00441449138559817,
            -0.0026916840713611884,
            -0.0021433533182417105,
            -0.011393931055161136,
            0.00831942078491083,
        ],
        [
            0.9924953570858336,
            -0.006527892583528865,
            0.00027703123670332134,
            -0.0043084140948845,
            -0.0025418940141105554,
            0.010403891427627623,
        ],
        [
            0.9910796745451873,
            -0.006066346677634581,
            0.0032124391214571626,
            -0.004684095647243339,
            0.005271523936204335,
            0.008559052737265285,
        ],
        [
            0.9898744236327127,
            -0.003589957116708843,
            0.005258145660315078,
            -0.002983415395102339,
            0.010601178931901838,
            0.004737087172996553,
        ],
        [
            0.9894058633306714,
            1.0720548026638724e-9,
            0.005986080494545013,
            -7.646033252523916e-9,
            0.012509735977239008,
            2.992840604048945e-9,
        ],
    ]

    stable = [converge(ic[1:6], reshape(ic[7:end], 6, 6), Φ; eps = 1e-8) for ic in ics]

    @test stable ≈ [
        [
            0.9894058633306991,
            -1.0720473066993234e-9,
            0.005986080494541589,
            7.645981012462348e-9,
            0.012509735977210676,
            -2.992820395381124e-9,
        ],
        [
            0.9898744236327279,
            0.003589957116710926,
            0.005258145660315757,
            0.002983415395071781,
            0.010601178931892018,
            -0.004737087172988469,
        ],
        [
            0.9910796745451955,
            0.006066346677635056,
            0.003212439121459165,
            0.004684095647225347,
            0.0052715239362023545,
            -0.008559052737263916,
        ],
        [
            0.9924953570858377,
            0.006527892583529262,
            0.0002770312367051028,
            0.004308414094874181,
            -0.0025418940141098013,
            -0.010403891427629895,
        ],
        [
            0.9934985344354048,
            0.004414491385598826,
            -0.0026916840713603093,
            0.00214335331823643,
            -0.011393931055160662,
            -0.008319420784914173,
        ],
        [
            0.9938017980401521,
            -2.637172844206485e-9,
            -0.004091057288364902,
            6.142460966514267e-9,
            -0.016385848085082626,
            5.354488437345373e-9,
        ],
        [
            0.9934985353164989,
            -0.00441449621708523,
            -0.002691681954481052,
            -0.0021433410126145015,
            -0.011393922984824547,
            0.008319427338101697,
        ],
        [
            0.9924953583827175,
            -0.006527896283858341,
            0.0002770340761640776,
            -0.004308400722454952,
            -0.002541885219410229,
            0.010403891704542608,
        ],
        [
            0.9910796756645397,
            -0.006066349416238057,
            0.003212441327022454,
            -0.004684081166596076,
            0.005271530075235392,
            0.008559049079820144,
        ],
        [
            0.9898744242561389,
            -0.00358995938694182,
            0.005258146784411463,
            -0.0029834002997269726,
            0.01060118194201734,
            0.004737081701461207,
        ],
        [
            0.9894058633306714,
            -1.0720398109166672e-9,
            0.005986080494545012,
            7.645928773607352e-9,
            0.012509735977239006,
            -2.9928001869866915e-9,
        ],
    ]
end

end