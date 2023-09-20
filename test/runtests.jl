"""
Restricted Two-body Model tests.
"""
module AstrodynamicalCalculationsTests

using Test, AstrodynamicalCalculations
using StaticArrays

@testset verbose = false "R2BP Determination" begin
    r = [0.0, 11681.0, 0.0]
    v = [5.134, 4.226, 2.787]
    μ = 398600.4354360959

    e, a, i, Ω, ω, ν = cartesian_to_keplerian(r, v, μ)
    @test isapprox(SVector(e, a, i, Ω, ω, ν),
        SVector(0.723452708202361,
            24509.265399338536,
            deg2rad(151.50460766373865),
            π / 2,
            deg2rad(270.0034742609256),
            deg2rad(89.99652573907436)),
        atol=1e-6)

    rₙ, vₙ = keplerian_to_cartesian(e, a, i, Ω, ω, ν, μ)
    @test isapprox(vcat(r, v), vcat(rₙ, vₙ), atol=1e-3)
end

@testset verbose = false "Kepler's Algorithm" begin
    r = [0.0, 11681.0, 0.0]
    v = [5.134, 4.226, 2.787]
    μ = 398600.4354360959
    a = semimajor_axis(r, v, μ)
    T = orbital_period(a, μ)

    rₙ, vₙ = kepler(r, v, μ, T)

    @test isapprox(vcat(r, v), vcat(rₙ, vₙ), atol=1e-3)
end

@testset verbose = true "Lambert Solvers" begin
    @testset "Universal" begin
        r = [0.0, 11681.0, 0.0]
        v = [5.134, 4.226, 2.787]
        Δt = 1000
        μ = 398600.4354360959

        rₙ, vₙ = kepler(r, v, μ, Δt; atol=1e-3)

        v₁, v₂ = lambert(r, rₙ, μ, Δt; trajectory=:short, atol=1e-6)

        @test isapprox(vcat(v, vₙ), vcat(v₁, v₂), atol=1e-3)
    end

    @testset "Lancaster / Blanchard" begin
        r = [0.0, 11681.0, 0.0]
        v = [5.134, 4.226, 2.787]
        Δt = 1000
        μ = 398600.4354360959

        rₙ, vₙ = kepler(r, v, μ, Δt; atol=1e-12)
        v₁, v₂ = R2BPCalculations.lambert_lancaster_blanchard(r,
            rₙ,
            Δt,
            μ;
            trajectory=:short,
            atol=1e-6)

        @test_broken isapprox(vcat(v, vₙ), vcat(v₁, v₂), atol=1e-3)
    end
end

@testset "CR3BP Calculations" begin
    r = [1.007988
        0.0
        0.001864]

    v = zeros(3, 1)

    μ = 3.003480593992993e-6

    @test jacobi_constant(r, v, μ) ≈ 3.000907212196274

    #
    # Conditions for tests below:
    # julia> μ = 0.012150584395829193
    # julia> u, T = halo(μ, 1; amplitude=0.005)
    #

    Φ = [1317.472125300217 -339.26725145920585 -22.23471304682866 388.2455372345198 126.49047353668445 -3.3779594177227; -429.7764385389933 111.53280315746214 7.269431052432993 -126.49047353678779 -41.404653982215095 1.0977750840404263; -11.440384647744368 2.9348175664641527 1.1929568568249131 -3.3779594177162653 -1.0977750840374354 -0.05845992955397675; 3612.12986893901 -929.3271447079647 -60.998106970610365 1064.4911782262984 346.9671305741014 -9.244834479682378; -1482.5514995781887 382.0700769975138 25.039782109090023 -437.2238230101127 -141.4481439160884 3.8211012689780377; -75.5369690753493 19.429643984545518 1.2797982252155324 -22.23471304678537 -7.269431052412963 1.1929568568243507]

    ics = [[0.823388563881332, 0.0, 0.005553604696093592, 0.0, 0.12683910108732768, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0], [0.8273403859561881, 0.03259254695615611, 0.004655244516514311, 0.027320231122920834, 0.10312408143143087, -0.006368351099157183, 1.3786797126365464, -0.08824308385714617, -0.009143474674588519, 2.826623191577671, -0.9767705092124925, -0.06678667433841168, -0.041844969630018845, 0.8820353502805051, 0.0007013669623268068, -0.45278244255732314, -0.7835124078992787, 0.0073360448800629365, -0.008855612000587742, 0.002327804586816214, 0.838573983007699, -0.06290404981028867, 0.02486880981909553, -1.144418551135245, 0.29476786414533784, -0.07907382422170534, -0.0008101636936140204, 1.2265481296819574, -0.604917980406487, -0.008685295114209327, 0.07388202008165795, 0.2492072545027064, -1.4738222079088099e-5, 0.5293937831553415, 0.7281497497436548, -0.00023580443544974562, -0.0007858587607819898, 0.00020348487224860056, 0.25939140376611963, -0.008269842352818602, 0.0028734095281463595, 0.8384737725375465], [0.8369339005610434, 0.05308879581928126, 0.002245931439623135, 0.03925760359490229, 0.04164563292173339, -0.010722629489038576, 2.633666680439342, -0.7476233640873686, -0.036115619838154166, 6.59317685013257, -4.305178122455261, -0.12567569616107588, -0.32236550652146106, 0.6479359223648997, 0.004637853579299471, -1.718449648616434, -0.7249471649546353, 0.021058591017896974, -0.03358790594797937, 0.01714406169353313, 0.40552901171028016, -0.1175501978131001, 0.08850970435530003, -1.9276108228616495, 0.7232939839971113, -0.36584119809468774, -0.005825038522167102, 2.0164150719091194, -1.588312972089223, -0.02821711404806066, 0.2835849868569326, 0.3549141154538107, -0.00031947991899352273, 0.9912017667556496, -0.033028457627487386, -0.0024977304687910985, -0.005384939024007069, 0.00272987184176848, 0.4345059352772165, -0.025760320739897816, 0.01772446341402129, 0.3998442189060235], [0.8470241424999021, 0.05395424966031193, -0.0009080043524526027, 0.031419954303289845, -0.036048701160724284, -0.011590792835804302, 5.2751979177709245, -2.8002759306940206, -0.06970307574708806, 13.575834570928198, -11.52199157988803, -0.08191686885532921, -1.0650447822856468, 0.6234811067801149, 0.010733752984669538, -3.905785161291004, 0.866145889611684, 0.016124697957608767, -0.07643671410242639, 0.05458515963145915, -0.16164597656191948, -0.2095060221160001, 0.19489307781221327, -2.0855174774020258, 1.5018686968935528, -1.0383366134481509, -0.014310928273230595, 3.9574099493625914, -3.5470463790411046, -0.02296768894075728, 0.6253995245350682, 0.18722436787118188, -0.0013253485604900028, 1.5577619859991105, -1.2832341322229566, -0.0027644056993001167, -0.014907405477056719, 0.010799843345176774, 0.4631273703125115, -0.043829837552022737, 0.04200385021401396, -0.20848785169720335], [0.8532353425901987, 0.03409350294370255, -0.003693332419778437, 0.0136102558436333, -0.10493897930691469, -0.007897889399969823, 11.061936482274414, -7.689747752640067, -0.03153665012283715, 31.785749672658245, -25.44043221945118, 0.5332758341624475, -2.6858333294059364, 1.3445895697563528, 0.002854429918083852, -8.652720643412609, 4.866690329415817, -0.113864600921434, -0.16640832220212004, 0.1333566308356175, -0.6638852788093453, -0.5058642813219177, 0.40239589734920744, -1.431763849780302, 3.195745363280542, -2.4844942161704124, -0.00386208695123994, 9.339460840937965, -7.368022339704923, 0.1499688501028994, 1.223461957527104, -0.41159173922619013, 0.002536978196004755, 3.113085020093486, -3.2037573059636526, 0.04850882843406646, -0.0306140680948119, 0.02606395260902634, 0.3190695726681844, -0.07676967973326097, 0.06982698735410356, -0.823253631643632], [0.854955140499413, 2.2843913370187256e-15, -0.004841260493981484, -8.426065623187225e-15, -0.1344033872050793, -2.3393929523281094e-16, 25.670001169891876, -17.290935930959833, 0.3735915323822458, 82.56024329676991, -44.66052062847337, 2.744550087867072, -6.549914989117929, 3.515781165459679, -0.09045012326731859, -21.560863766695647, 11.122750444783911, -0.65681735247777, -0.4081906725844075, 0.28639114669384175, -0.8777886212624451, -1.3960821106365155, 0.7204334759076325, -0.04847727067031394, 7.495647097356273, -5.235090649386511, 0.11199867561047291, 24.307495565860556, -12.730790263231905, 0.7894471966986759, 2.6328847934916584, -1.602867310555306, 0.041320768129938484, 7.950197124746471, -5.4159521324395765, 0.26854701213925986, -0.06452008851439175, 0.049418531406256176, 0.03736159910558069, -0.19318325888690766, 0.10035816116241962, -1.1539010608138158], [0.8532353425901938, -0.03409350294369572, -0.0036933324197787306, -0.013610255843655042, -0.10493897930691656, 0.007897889399968805, 61.783725674727414, -31.809460366393743, 1.5639194613448177, 192.5622814708414, -60.811745821785316, 5.790230207958852, -15.9467057502335, 7.366098751441629, -0.3866391287246332, -50.036359202730644, 16.771181856064608, -1.474065945922604, -1.0242029826389052, 0.5249458158655891, -0.6919274956044533, -3.288510633984684, 1.0165171629160492, 1.321471791772166, 18.12710682888481, -9.376464402158039, 0.4560226259807858, 56.69804470264641, -17.397205275984405, 1.6789260061002176, 6.108694482836869, -3.292222472289838, 0.15851822297241477, 18.49715015775276, -6.775733307775279, 0.5687877877351285, -0.1529731671956582, 0.08076178851334609, -0.2644398496513836, -0.4861071061444217, 0.12936352157783781, -0.9612148163601355], [0.8470241424998886, -0.053954249660305716, -0.0009080043524529912, -0.031419954303331776, -0.036048701160720203, 0.011590792835803694, 139.81168166322985, -51.900248046349844, 3.271613747195812, 397.8858051004851, -91.85758770599232, 5.821945264799484, -36.17535833489009, 12.967888420563526, -0.8279061449068151, -102.92655322720218, 25.534744119150396, -1.5228161048932043, -2.3513343309735673, 0.866235468830786, -0.2198726250639241, -6.7447477075950895, 1.5755137371116215, 1.987134962713231, 41.11521926552655, -15.166743448559673, 0.9526780029549172, 117.27289978925529, -26.666892209713573, 1.6975459397425652, 13.579339684552384, -5.4115134792664845, 0.32466549527180183, 38.02628330961531, -9.24871447769774, 0.5601717489178005, -0.35373005498093574, 0.12514521971200396, -0.45673669921893095, -1.0315564554873287, 0.2132601376815499, -0.4067512272598734], [0.8369339005610131, -0.05308879581927243, 0.0022459314396227116, -0.03925760359498695, 0.04164563292175057, 0.010722629489038737, 298.1536849578567, -88.79218713147596, 4.108212352069355, 805.8842647516645, -197.03199719416259, -1.3587839120057053, -77.0304035805354, 23.008804194914685, -1.051997832312519, -207.58784021007904, 52.69517331340727, 0.31433070459529, -5.028116277146918, 1.4986245228726873, 0.3336885466879633, -13.604470125477365, 3.3661481657079335, 1.9537305723180947, 87.79671835004706, -25.956894579473673, 1.196556663820022, 237.60479660840883, -57.88401149438683, -0.4014318236821094, 28.712231594102008, -8.994281477261307, 0.40346490381046596, 77.0693721005221, -18.784314465348764, -0.14157890224819605, -0.7638066834236382, 0.21608393326675798, -0.4802205363194697, -2.081496212262708, 0.5023315171518498, 0.23330573472502658], [0.8273403859561236, -0.032592546956138974, 0.004655244516514174, -0.027320231123097893, 0.1031240814314825, 0.006368351099159685, 624.1563531688901, -178.97128138129779, 1.2034561562153203, 1685.807676266411, -518.1477491379977, -23.313639430926575, -160.92062399834495, 46.705824906057714, -0.31350301921564855, -433.6254768006774, 134.82826866054083, 5.971159098990085, -10.529645901486171, 3.0305772015792436, 0.8168304281567269, -28.450876102792893, 8.77375084473596, 1.5423565889576247, 183.90492595796482, -52.51978660929593, 0.34248355473230047, 496.93345370466517, -152.78832660830577, -6.855854030546668, 59.93604732520242, -17.543497122306604, 0.12072458185754525, 161.65110999609013, -49.108411393866255, -2.2542758104173393, -1.6023417655780476, 0.4493519667693413, -0.33589553197288335, -4.324479481929502, 1.3446971876129539, 0.797006208427492], [0.8233885638811947, 4.33547133764043e-14, 0.00555360469609479, -3.778286963052016e-13, 0.12683910108748223, 7.86884919269592e-15, 1317.472125300217, -429.7764385389933, -11.440384647744368, 3612.12986893901, -1482.5514995781887, -75.5369690753493, -339.26725145920585, 111.53280315746214, 2.9348175664641527, -929.3271447079647, 382.0700769975138, 19.429643984545518, -22.23471304682866, 7.269431052432993, 1.1929568568249131, -60.998106970610365, 25.039782109090023, 1.2797982252155324, 388.2455372345198, -126.49047353678779, -3.3779594177162653, 1064.4911782262984, -437.2238230101127, -22.23471304678537, 126.49047353668445, -41.404653982215095, -1.0977750840374354, 346.9671305741014, -141.4481439160884, -7.269431052412963, -3.3779594177227, 1.0977750840404263, -0.05845992955397675, -9.244834479682378, 3.8211012689780377, 1.1929568568243507]]

    unstable = [diverge(ic[1:6], reshape(ic[7:end], 6, 6), Φ; eps=1e-8) for ic in ics]

    @test unstable ≈ [[0.8233885670604992, -1.0366541345478615e-9, 0.005553604668461628, 8.715859137046784e-9, 0.12683909750709388, -1.8222484329190235e-10], [0.8273403889974105, 0.03259254570178936, 0.004655244468062643, 0.027320239369938236, 0.1031240768350659, -0.006368351287082361], [0.8369339034694137, 0.05308879425915669, 0.0022459313856777996, 0.039257611236035383, 0.04164562738118049, -0.010722629629902078], [0.8470241452302693, 0.05395424782876832, -0.0009080043935676615, 0.031419961538898866, -0.036048707229790644, -0.011590792871138664], [0.8532353451489048, 0.03409350102660894, -0.0036933324275151497, 0.013610263353121406, -0.10493898508343963, -0.00789788926934211], [0.8549551430909017, -1.7690279810307236e-9, -0.0048412604557870696, 8.408547056141493e-9, -0.13440339160640008, 2.7897965198713934e-10], [0.8532353454774562, -0.03409350442176633, -0.003693332346613867, -0.013610246811438523, -0.104938982103295, 0.007897889670159345], [0.8470241457112383, -0.053954250841694876, -0.000908004277468657, -0.03141994514599701, -0.036048703261041745, 0.01159079296921633], [0.8369339039236694, -0.05308879681438595, 0.0022459314858005795, -0.03925759449829196, 0.041645630696876774, 0.010722629473653931], [0.8273403892761712, -0.032592547905735625, 0.004655244522837673, -0.0273202221542787, 0.10312407867097084, 0.0063683509751929325], [0.823388567060362, -1.036610779834484e-9, 0.005553604668462826, 8.71548130835048e-9, 0.12683909750724842, -1.8221697444270972e-10]]

    stable = [converge(ic[1:6], reshape(ic[7:end], 6, 6), Φ; eps=1e-8) for ic in ics]

    @test stable ≈ [[0.8233885607021647, -1.036654134541258e-9, 0.005553604723725555, 8.715859137047568e-9, 0.1268391046675615, -1.8222484329220807e-10], [0.8273403826361405, 0.03259254600655946, 0.004655244510190812, 0.027320240091740028, 0.10312408419194252, -0.0063683512231239355], [0.8369338971983872, 0.05308879482416774, 0.002245931393445267, 0.03925761269159728, 0.041645635146607185, -0.010722629504423382], [0.8470241392885524, 0.05395424847892277, -0.0009080044274369369, 0.03141996346062461, -0.03604869906040274, -0.011590792702391666], [0.8532353397029364, 0.03409350146563194, -0.0036933324929433005, 0.013610264875849818, -0.10493897651053626, -0.007897889129779283], [0.8549551379079243, -1.7690279812603002e-9, -0.004841260532175899, 8.408547056214322e-9, -0.13440338280375852, 2.7897965200664483e-10], [0.8532353400314878, -0.034093504860789325, -0.003693332412042018, -0.013610248334166937, -0.10493897353039162, 0.007897889530596518], [0.8470241397695214, -0.05395425149184933, -0.0009080043113379323, -0.03141994706772276, -0.03604869509165384, 0.011590792800469332], [0.8369338976526429, -0.053088797379397, 0.0022459314935680467, -0.039257595953853856, 0.04164563846230347, 0.010722629348175236], [0.8273403829149012, -0.03259254821050572, 0.004655244564965842, -0.027320222876080488, 0.10312408602784746, 0.006368350911234508], [0.8233885607020275, -1.0366107759879842e-9, 0.005553604723726753, 8.715481299812045e-9, 0.12683910466771603, -1.8221697411827795e-10]]
end

end # module