using MassJ, Test
using Plots
using DataStructures

function tests()
    @testset "Subset of tests"  begin
        inf = MassJ.info("test.mzXML", verbose = true)
        @test inf[1] == "parentFile: test.raw"                                         #1
        @test inf[9] == "6 scans"                                                      #2
        @test inf[10] == "MS1+"                                                        #3
        @test inf[11] == "MS2+ 1255.5  CID(CE=18)"                                     #4
        @test inf[12] == "MS3+ 902.33  PQD(CE=35)"                                     #5

        scans = MassJ.load("test.mzXML")
        @test eltype(scans)              == MassJ.MSscan                                 #6
        @test length(scans)              == 6                                          #7
        @test scans[1].num               == 1                                          #8
        @test scans[2].level             == 2                                          #9
        @test scans[3].polarity          == "+"                                        #10
        @test scans[2].activationMethod  == "CID"                                      #11
        @test scans[3].collisionEnergy   == 35.0                                       #12
        @test size(scans[1].int, 1)      == 22320                                      #13

        rt = MassJ.retention_time("test.mzXML")
        @test length(rt) == 6                                                          #14

        cr = MassJ.chromatogram("test.mzXML", method = MassJ.TIC() )
        @test length(cr.rt) == 6                                                       #15

        cr = MassJ.chromatogram("test.mzXML", method = MassJ.MZ([0, 500]))
        @test length(cr.rt) == 6                                                       #16

        cr = MassJ.chromatogram("test.mzXML", method = MassJ.∆MZ([1000, 1]))
        @test length(cr.rt) == 6                                                       #17

        cr = MassJ.chromatogram("test.mzXML", method = MassJ.BasePeak())
        @test length(cr.rt) == 6                                                       #18

        cr = MassJ.chromatogram("test.mzXML", MassJ.Polarity("+"), MassJ.Scan(2),MassJ.Precursor(1255.5), MassJ.Activation_Energy(18), MassJ.Activation_Method("CID"), MassJ.Level(2) )
        @test length(cr.rt) == 1                                                       #19

        rt = MassJ.retention_time(scans)
        @test length(rt) == 6                                                          #20

        cr = MassJ.chromatogram(scans, method = MassJ.TIC() )
        @test length(cr.rt) == 6                                                       #21

        cr = MassJ.chromatogram(scans, method = MassJ.MZ([0, 500]))
        @test length(cr.rt) == 6                                                       #22

        cr = MassJ.chromatogram(scans, method = MassJ.∆MZ([1000, 1]))
        @test length(cr.rt) == 6                                                       #23

        cr = MassJ.chromatogram(scans, method = MassJ.BasePeak())
        @test length(cr.rt) == 6                                                       #24

        cr = MassJ.chromatogram(scans, MassJ.Polarity("+"),MassJ.Scan(2),MassJ.Precursor(1255.5),MassJ.Activation_Energy(18),MassJ.Activation_Method("CID"),MassJ.Level(2) )
        @test (cr.rt, cr.ic) == ([0.7307], [9727.2])                                   #25

        cr = MassJ.chromatogram(scans, MassJ.Polarity(["+"]),MassJ.Scan([2,3]),MassJ.Precursor([1255.5, 902.33]),MassJ.Activation_Energy([18, 35]),MassJ.Activation_Method(["CID", "PQD"]),MassJ.Level([2, 3]) )
        @test (cr.rt, cr.ic) == ([0.7307, 2.1379], [9727.2, 11.3032])                  #26

        ms = MassJ.average("test.mzXML")
        @test length(ms.num) == 6                                                      #27

        ms = MassJ.average("test.mzXML", MassJ.Polarity("+"),MassJ.Scan(2),MassJ.Precursor(1255.5),MassJ.Activation_Energy(18),MassJ.Activation_Method("CID"),MassJ.RT(1),MassJ.IC([0, 1e4]))
        @test ms isa MassJ.MSscan                                                        #28
        @test ms.num == 2                                                              #29

        ms = MassJ.average(scans)
        @test length(ms.num) == 6                                                      #30

        ms = MassJ.average(scans, MassJ.Polarity("+"),MassJ.Scan(2),MassJ.Precursor(1255.5),MassJ.Activation_Energy(18),MassJ.Activation_Method("CID"),MassJ.RT(1),MassJ.IC([0, 1e4]))
        @test ms isa MassJ.MSscan                                                        #31
        @test ms.num == 2                                                              #32

        ms = MassJ.average(scans, MassJ.Polarity(["+"]),MassJ.Scan([2,3]),MassJ.Precursor([1255.5, 902.33]),MassJ.Activation_Energy([18, 35]),MassJ.Activation_Method(["CID", "PQD"]),MassJ.RT([1,2]),MassJ.IC([0, 1e4]))
        @test ms isa MassJ.MSscans                                                       #33
        @test ms.num == [2, 3]                                                         #34

        ms = MassJ.average("test.mzXML", MassJ.RT( [[1,2], [2,3]] ), stats = false )
        @test ms isa MassJ.MSscans                                                       #35
        @test ms.num == [2, 3, 4]                                                      #36

        ms = MassJ.average("test.mzXML", MassJ.Polarity(["+"]),MassJ.Scan([2,3]),MassJ.Precursor([1255.5, 902.33]),MassJ.Activation_Method(["CID", "PQD"]),MassJ.RT([1,2]),MassJ.IC([0, 1e4]))   #MassJ.Activation_Energy([18., 35.]),
        @test ms isa MassJ.MSscans                                                       #37
        @test ms.num == [2, 3]                                                         #38

        cr = MassJ.chromatogram("test.mzXML", MassJ.Polarity(["+"]),MassJ.Scan([2,3]),MassJ.Precursor([1255.5, 902.33]),MassJ.Activation_Method(["CID", "PQD"]),MassJ.Level([2, 3]) )   #MassJ.Activation_Energy([18.0, 35.0]),
        @test length(cr.rt) == 2                                                       #39

        ms = MassJ.average(scans, MassJ.RT( [[1,2], [2,3]] ), stats = false )
        @test ms isa MassJ.MSscans                                                       #40
        @test ms.num == [2, 3, 4]                                                      #41

        a = scans[1] / 2.
        @test a.tic == 2.540975e6                                                      #42

        a = scans[1] * 2.
        @test a.tic == 1.01639e7                                                       #43

        a = ms * 2.
        @test a.tic == 3.2120923354666666e6                                            #44

        a = 2. * scans[1]
        @test a.tic == 1.01639e7                                                       #45

        a = scans[1] * scans[2]
        @test a.tic == 4.943314404e10                                                  #46

        a = scans[2] * scans[1]
        @test a.tic == 4.943314404e10                                                  #47

        a = scans[1] - scans[2]
        @test a.tic == 5.0722228e6                                                     #48

        a = scans[2] - scans[1]
        @test a.tic == -5.0722228e6                                                    #49

        a = scans[1] - scans[4]
        @test a.tic == 273550.0                                                        #50

        b = ms - scans[1]
        @test b.num == [2,3,4]                                                         #51

        b = ms + scans[1]
        @test b.num == [2,3,4,1]                                                       #52

        b = scans[1] + ms
        @test b.num == [1,2,3,4]                                                       #53

        b = scans[1] + scans[4]
        @test b.num == [1,4]                                                           #54

        a = MassJ.avg(scans[1], scans[2])
        @test a.num == [1,2]                                                           #55

        a = MassJ.avg(scans[1], scans[4])
        @test a.num == [1,4]                                                           #56

        info = MassJ.info("test64.mzXML")
        @test info[2] == "MS1-"                                                        #57

        scans = MassJ.load("test64.mzXML")
        @test eltype(scans)              == MassJ.MSscan                                 #58

        info = MassJ.info("test.mzXMLM")
        @test info.msg == "File format not supported."                                 #59

        scans = MassJ.load("test.mzXMLL")
        @test info.msg == "File format not supported."                                 #60

        scans = MassJ.load("bad1.mzXML")
        @test scans.msg == "Not an mzXML file."                                        #61

        scans = MassJ.info("bad1.mzXML")
        @test scans.msg == "Not an mzXML file."                                        #62

        scans = MassJ.load("bad2.mzXML")
        @test scans[1].num == 0                                                        #63

        scans = MassJ.load("bad3.mzXML")
        @test scans[1].num == scans[2].num == 0                                        #64

        cr = MassJ.chromatogram("test.mzXML", method = MassJ.∆MZ([1, 2]))
        @test cr.msg == "Bad mz ± ∆mz values."                                         #65

        cr = MassJ.chromatogram(scans, method = MassJ.∆MZ([1, 2]))
        @test cr.msg == "Bad mz ± ∆mz values."                                         #66

        scans = MassJ.load("test.mzXML")
        @test MassJ.smooth(scans[1], method = MassJ.SG(7,15,0)) isa MassJ.MSscan             #67

        a = MassJ.avg(scans[1], scans[4])
        @test MassJ.smooth(a) isa MassJ.MSscans                                            #68

        a = scans[1] * scans[4]
        @test a.num == [1,4]                                                           #69

        a = (scans[2]+scans[3]) - scans[1]
        @test a.num == [2, 3]                                                          #70

        a = (scans[1] + scans[4]) - (scans[1] - scans[4])
        @test a.num == [1, 4]                                                          #71

        a = scans[1] + MassJ.avg(scans[2], scans[5])
        @test a.num == [1, 2, 5]                                                       #72

        a = MassJ.smooth(scans[1], method = MassJ.SG(5,9,0))
        @test a.num == 1                                                               #73

       a = MassJ.centroid(scans[1], method = MassJ.TBPD(:gauss, 4500., 0.2))               #74
       @test length(a.int) == 957

       @test typeof(plot(scans[1], method = :relative)) == Plots.Plot{Plots.GRBackend} #75
       @test typeof(plot(scans[1], method = :absolute)) == Plots.Plot{Plots.GRBackend} #76

       a = MassJ.avg(scans[2], scans[5])
       @test typeof(plot( a, method = :relative )) == Plots.Plot{Plots.GRBackend}      #77
       @test typeof(plot( a, method = :absolute )) == Plots.Plot{Plots.GRBackend}      #78

       cr = MassJ.chromatogram(scans)
       @test typeof(plot( cr, method = :relative )) == Plots.Plot{Plots.GRBackend}     #79
       @test typeof(plot( cr, method = :absolute )) == Plots.Plot{Plots.GRBackend}     #80

       a = MassJ.centroid(scans[1], method = MassJ.TBPD(:voigt, 4500., 0.2))
       @test length(a.int) == 961                                                      #81

       a = MassJ.centroid(scans[1], method = MassJ.TBPD(:lorentz, 4500., 0.2))
       @test length(a.int) == 964                                                      #82

       a = MassJ.centroid(scans[1], method = MassJ.TBPD(:other, 4500., 0.2))
       @test a.msg == "Unsupported peak profile. Use :gauss, :lorentz or :voigt."      #83

       a = MassJ.centroid(scans[1], method = MassJ.SNRA(1., 100))
       @test length(a.int) == 109                                                      #84

       s1 = MassJ.extract(scans, MassJ.Activation_Energy([18,35]))
       @test length(s1) == 4                                                           #85

       s1 = MassJ.extract("test.mzXML", MassJ.Activation_Energy(18))
       @test length(s1) == 2                                                           #86

       s1 = MassJ.extract(scans, MassJ.Scan(1))
       @test length(s1) == 1                                                           #87

       s1 = MassJ.extract("test.mzXML", MassJ.Scan(1))
       @test length(s1) == 1                                                           #88

       bs = MassJ.baseline_correction(scans, method = MassJ.TopHat(1))
       @test length(bs) == 6                                                           #89

       bs = MassJ.baseline_correction(scans[1], method = MassJ.TopHat(1))
       @test length(bs.int) == length(scans[1].int)                                    #90

       c = MassJ.centroid(scans, method = MassJ.TBPD(:gauss, 4500., 0.2)) ;
       bs = MassJ.baseline_correction( c, method = MassJ.LOESS(3))
       @test length(bs) == 6                                                           #91

       bs = MassJ.baseline_correction(c[1], method = MassJ.LOESS(3))
       @test length(bs.int) == length(c[1].int)                                        #92

       bs = MassJ.baseline_correction(scans, method = MassJ.IPSA(51,100))
       @test length(bs) == 6                                                           #93

       bs = MassJ.baseline_correction(scans[1], method = MassJ.IPSA(51,100))
       @test length(bs.int) == length(scans[1].int)                                    #94

       a = MassJ.smooth(scans, method = MassJ.SG(5,9,0))
       @test length(a) == 6                                                            #95

       c = MassJ.centroid(scans, method = MassJ.TBPD(:lorentz, 4500., 0.2)) ;
       d = MassJ.centroid(scans, method = MassJ.TBPD(:voigt, 4500., 0.2)) ;
       @test length(c) == length(d)                                                    #96

       a = MassJ.centroid(scans[3], method = MassJ.SNRA(1., 100))
       @test length(a.int) == 0                                                        #97

       bs = MassJ.baseline_correction(scans[1], method = MassJ.IPSA(50,100))
       @test length(bs.int) == length(scans[1].int)                                    #98

       cr = MassJ.chromatogram(scans, method = MassJ.BasePeak() )
       @test length(cr.rt) == 6                                                        #99

       f = MassJ.formula("CH3(13C)10H3Kr(NaH2)2")                                        #100
      @test f == Dict("Na" => 2,"Kr" => 1,"C" => 1,"13C" => 10,"H" => 10)

       m = MassJ.masses("C254 H377 N65 O75 S6")                                          #101
      @test m["Monoisotopic"] ≈ 5729.60087099839
      @test m["Average"] ≈ 5733.55
      @test m["Nominal"] ≈ 5727.0

      I = MassJ.isotopic_distribution("CH4", 0.9999, charge = +1)                        #102
      @test I[2,1:end] == [16.03130012908, 0.9887541751052761, 1, 0, 4, 0]

      a = MassJ.simulate(I, 0.4, Npoints = 5)                                            #103
      @test a.int == [100.0, 6.035851011021856, 0.06994625998243831, 1.1368602493290418, 0.06861901924947518]

       m = MassJ.masses(f)                                                               #104
      @test m == Dict("Monoisotopic" => 282.0028349717, "Average" => 281.902086912, "Nominal" => 282.0)

      I = MassJ.isotopic_distribution(f, 0.9999, charge = +1)                            #105
      @test I[2,1:end][1:2] ≈ [282.0028349717, 0.5630635281692917]

      a = MassJ.simulate(I, 0.4, model=:lorentz, Npoints = 5)                            #106
      @test a.int ≈ [0.5359868694750152, 16.332387108915455, 100.0, 0.560482304632663, 0.0834376623225204]

    end
end


function test_isotopes()
    @testset "Isotopes - type stability and optimizations" begin

        # Elements dict is const (enables type inference)
        @test isconst(MassJ, :Elements)

        # PriorityQueues are typed
        pq = PriorityQueue{Vector{Int},Float64}()
        pq[[1,0]] = 0.5
        pq[[0,1]] = 0.3
        @test peek(pq) == ([0,1] => 0.3)   # min-heap: lowest first

        # formula parsing
        @test MassJ.formula("C2H6O") == Dict("C" => 2, "H" => 6, "O" => 1)
        @test MassJ.formula("H2O") == Dict("H" => 2, "O" => 1)
        @test MassJ.formula("NaCl") == Dict("Na" => 1, "Cl" => 1)
        @test_throws ErrorException MassJ.formula("123bad")
        @test_throws ErrorException MassJ.formula("Xx")

        # masses calculation
        m = MassJ.masses("H2O")
        @test m["Monoisotopic"] ≈ 18.01056468474
        @test m["Average"] ≈ 18.015
        @test m["Nominal"] == 18.0

        m2 = MassJ.masses(Dict("H" => 2, "O" => 1))
        @test m2 == m

        # stirling approximation vs exact log factorial
        @test MassJ.stirling(100) ≈ log(factorial(big(100))) atol=0.01
        @test MassJ.stirling(500) > 0

        # isotopologue_probability - low mass path (Natoms < 20)
        prob_H2 = MassJ.isotopologue_probability(Dict("H" => 2), Dict("H" => [2, 0]), MassJ.Elements)
        @test prob_H2 ≈ MassJ.Elements["H"][1].f^2 atol=1e-10

        # isotopologue_probability - high mass path (Natoms >= 20)
        prob_C20 = MassJ.isotopologue_probability(Dict("C" => 20), Dict("C" => [20, 0]), MassJ.Elements)
        @test prob_C20 ≈ 0.8049835604738165 atol=1e-8

        # isotopologue_mass
        cm = MassJ.isotopologue_mass([Pair("H", 2)], Dict("H" => [2, 0]), MassJ.Elements)
        @test cm ≈ 2 * MassJ.Elements["H"][1].m

        # most_probable_isotopologue
        alpha = MassJ.most_probable_isotopologue(Dict("C" => 10, "H" => 22), MassJ.Elements)
        @test sum(alpha["C"]) == 10     # conservation
        @test sum(alpha["H"]) == 22     # conservation
        @test alpha["C"][1] >= alpha["C"][2]  # 12C more abundant than 13C

        # hill_climbing finds optimum
        f_obj(x) = -(x[1] - 3)^2 - (x[2] - 2)^2
        P = MassJ.hill_climbing([1, 4], f_obj)
        @test P == [3, 2]

        # hill_climbing with single-element vector (early return)
        P_single = MassJ.hill_climbing([5], x -> -x[1]^2)
        @test P_single == [5]  # no neighbors possible

        # isotopic distribution - basic
        I = MassJ.isotopic_distribution("H2O", 0.99, charge = 1)
        @test size(I, 2) == 7   # Masses, Probability, + isotope columns
        @test I[2, 1] ≈ 18.01056468474   # monoisotopic mass
        @test I[2, 2] ≈ 0.9973367663173334  # probability
        @test I[1, 1] == "Masses"
        @test I[1, 2] == "Probability"

        # isotopic distribution - charge state divides mass
        I2 = MassJ.isotopic_distribution("H2O", 0.99, charge = 2)
        @test I2[2, 1] ≈ 18.01056468474 / 2

        # isotopic distribution - larger molecule
        I3 = MassJ.isotopic_distribution("C254 H377 N65 O75 S6", 0.5)
        @test size(I3, 1) > 2          # more than just header
        @test I3[1, 1] == "Masses"
        probs = [I3[i, 2] for i in 2:size(I3, 1)]
        @test all(p -> p > 0, probs)   # all probabilities positive
        @test sum(probs) ≈ 0.5 atol=0.1

    end
end


function test_deconvolution()
    @testset "Deconvolution - helpers and integration" begin

        # is_evenly_spaced
        @test MassJ.is_evenly_spaced([1.0, 2.0, 3.0, 4.0]) == true
        @test MassJ.is_evenly_spaced([1.0, 2.0, 3.5, 4.0]) == false
        @test MassJ.is_evenly_spaced(collect(range(100.0, stop=200.0, length=100))) == true

        # resampling
        X = [1.0, 2.0, 3.5, 5.0]
        Y = [10.0, 20.0, 35.0, 50.0]
        newX, newY = MassJ.resampling(X, Y)
        @test length(newX) == length(X)
        @test newX[1] ≈ X[1]
        @test newX[end] ≈ X[end]
        @test MassJ.is_evenly_spaced(newX)

        # new_mass
        @test MassJ.new_mass(500.0, 5, 1, 1.00782503227) ≈ 416.83463750537834
        @test MassJ.new_mass(500.0, 5, -1, 1.00782503227) ≈ 624.7480437419325
        # new_mass with zero shift returns original
        @test MassJ.new_mass(500.0, 5, 0, 1.0) ≈ 500.0

        # get_peak_shape
        mz = collect(range(100.0, stop=200.0, length=1000))
        g = MassJ.get_peak_shape(MassJ.gauss, 0.5, mz)
        @test length(g) > 0
        @test sum(g) ≈ 1.0 atol=1e-10      # normalized
        @test g[argmax(g)] > 0              # has a peak
        # approximately symmetric shape
        @test g[1] ≈ g[end] rtol=1.0

        # get_peak_shape with lorentz
        g_lor = MassJ.get_peak_shape(MassJ.lorentz, 0.5, mz)
        @test sum(g_lor) ≈ 1.0 atol=1e-10

        # figure_of_merit
        h = [1.0, 2.0, 3.0, 4.0, 5.0]
        @test MassJ.figure_of_merit(h, h) ≈ 1.0           # perfect match
        @test MassJ.figure_of_merit(h, h .+ 0.1) ≈ 0.995  # close match
        @test MassJ.figure_of_merit(h, h .+ 0.1) > 0.99

        # chargefilter - basic shape and properties
        mz_test = collect(range(400.0, stop=600.0, length=200))
        int_test = zeros(200)
        int_test[100] = 100.0   # single peak at ~500
        f_test = zeros(10, 200)
        for i in 1:10
            f_test[i, :] = int_test
        end
        s = MassJ.chargefilter(mz_test, int_test, f_test, (1, 10), 1.0, 1)
        @test size(s) == (10, 200)
        @test all(s .>= 0)     # non-negative output

        # project_N_convolve
        g_small = MassJ.get_peak_shape(MassJ.gauss, 2.0, mz_test)
        c = MassJ.project_N_convolve(s, g_small)
        @test length(c) == 200
        @test all(c .>= 0)

        # Charges struct construction
        ch = MassJ.Charges(adduct="H", range=(5, 15), width=2)
        @test ch.adduct == "H"
        @test ch.range == (5, 15)
        @test ch.width == 2

        # Charges with default width
        ch2 = MassJ.Charges(adduct="Na", range=(1, 5))
        @test ch2.width == 1

        # _resolve_shape
        scans = MassJ.load("test.mzXML")
        model_gauss = MassJ._resolve_shape(scans[1], :gauss)
        @test model_gauss === MassJ.gauss
        model_lor = MassJ._resolve_shape(scans[1], :lorentz)
        @test model_lor === MassJ.lorentz
        model_voigt = MassJ._resolve_shape(scans[1], :voigt)
        @test model_voigt === MassJ.voigt
        @test_throws ErrorException MassJ._resolve_shape(scans[1], :invalid)

        # _resolve_FWHM with explicit values
        @test MassJ._resolve_FWHM(scans[1], MassJ.gauss, -1, 0.5) == 0.5   # explicit FWHM
        @test MassJ._resolve_FWHM(scans[1], MassJ.gauss, 1000.0, -1) == 0.5  # R=1000 → 500/1000=0.5
        @test_throws ErrorException MassJ._resolve_FWHM(scans[1], MassJ.gauss, 1000.0, 0.5)  # both

        # roughguess_FWHM runs without error
        w = MassJ.roughguess_FWHM(scans[1])
        @test w > 0

        # deconv method exists with correct dispatch
        @test hasmethod(MassJ.deconv, Tuple{MassJ.MSscan, MassJ.Charges})
        @test hasmethod(MassJ.deconv, Tuple{MassJ.MSscans, MassJ.Charges})

    end
end


function test_interpolation_import()
    @testset "Interpolations import fix" begin
        # Line is accessible from MassJ (needed for extrapolation)
        @test isdefined(MassJ, :Line)
        @test isdefined(MassJ, :LinearInterpolation)

        # Verify interpolation with Line extrapolation works
        x = [1.0, 2.0, 3.0, 4.0]
        y = [10.0, 20.0, 30.0, 40.0]
        itp = MassJ.LinearInterpolation(x, y, extrapolation_bc=MassJ.Line())
        @test itp(2.5) ≈ 25.0
        @test itp(0.0) ≈ 0.0    # extrapolated
    end
end


function test_mzml()
    @testset "mzML format" begin
        # info
        inf = MassJ.info("test.mzML")
        @test any(contains(s, "3 scans") for s in inf)
        @test any(contains(s, "MS1+") for s in inf)
        @test any(contains(s, "MS2+") for s in inf)

        inf_v = MassJ.info("test.mzML", verbose=true)
        @test any(contains(s, "test.raw") for s in inf_v)
        @test any(contains(s, "LTQ FT") for s in inf_v)

        # load all
        scans = MassJ.load("test.mzML")
        @test length(scans) == 3
        @test eltype(scans) == MassJ.MSscan

        # Scan 1: MS1+, profile, rt=0.5 min
        s1 = scans[1]
        @test s1.level == 1
        @test s1.rt ≈ 0.5
        @test s1.polarity == "+"
        @test s1.spectrumType == :profile
        @test s1.tic ≈ 19000.0
        @test s1.basePeakMz ≈ 400.0
        @test s1.basePeakIntensity ≈ 8000.0
        @test s1.precursor ≈ 0.0
        @test s1.chargeState == 0
        @test length(s1.mz) == 5
        @test s1.mz ≈ [100.0, 200.0, 300.0, 400.0, 500.0]
        @test s1.int ≈ [1000.0, 5000.0, 3000.0, 8000.0, 2000.0]

        # Scan 2: MS2+ CID, centroid, rt=1.0 min (60s converted)
        s2 = scans[2]
        @test s2.level == 2
        @test s2.rt ≈ 1.0
        @test s2.polarity == "+"
        @test s2.spectrumType == :centroid
        @test s2.precursor ≈ 400.0
        @test s2.chargeState == 2
        @test s2.activationMethod == "CID"
        @test s2.collisionEnergy ≈ 25.0
        @test length(s2.mz) == 4
        @test s2.mz ≈ [110.0, 150.0, 200.0, 250.0]

        # Scan 3: MS2+ HCD
        s3 = scans[3]
        @test s3.level == 2
        @test s3.rt ≈ 1.5
        @test s3.precursor ≈ 500.0
        @test s3.chargeState == 3
        @test s3.activationMethod == "HCD"
        @test s3.collisionEnergy ≈ 30.0
        @test length(s3.mz) == 3

        # retention_time
        rt = MassJ.retention_time("test.mzML")
        @test length(rt) == 3
        @test rt ≈ [0.5, 1.0, 1.5]

        # chromatogram (TIC)
        chrom = MassJ.chromatogram("test.mzML")
        @test length(chrom.rt) == 3
        @test chrom.rt ≈ [0.5, 1.0, 1.5]
        @test chrom.ic ≈ [19000.0, 4800.0, 2100.0]
        @test chrom.maxic ≈ 19000.0

        # chromatogram with filter
        chrom2 = MassJ.chromatogram("test.mzML", MassJ.Level(2))
        @test length(chrom2.rt) == 2

        # chromatogram base peak
        chrom_bp = MassJ.chromatogram("test.mzML", method=MassJ.BasePeak())
        @test chrom_bp.ic ≈ [8000.0, 2000.0, 1200.0]

        # extract
        ms2 = MassJ.extract("test.mzML", MassJ.Level(2))
        @test length(ms2) == 2
        @test ms2[1].precursor ≈ 400.0
        @test ms2[2].precursor ≈ 500.0

        # average
        avg = MassJ.average("test.mzML", MassJ.Level(2))
        @test avg isa MassJ.MSscans

        # New fields present
        @test s1.mobilityType == :none
        @test s1.driftTime ≈ -1.0
        @test s1.compensationVoltage ≈ 0.0
        @test s1.metadata isa Dict{String,Any}
    end
end


function test_mgf()
    @testset "MGF format" begin
        # info
        inf = MassJ.info("test.mgf")
        @test inf[1] == "3 scans"
        @test any(contains(s, "400.0") for s in inf)
        @test any(contains(s, "500.0") for s in inf)
        @test any(contains(s, "600.0") for s in inf)

        # load all
        scans = MassJ.load("test.mgf")
        @test length(scans) == 3
        @test eltype(scans) == MassJ.MSscan

        # Scan 1
        s1 = scans[1]
        @test s1.num == 1        # sequential index
        @test s1.rt ≈ 0.5       # 30s / 60
        @test s1.level == 2      # MGF default
        @test s1.precursor ≈ 400.0
        @test s1.chargeState == 2
        @test s1.polarity == "+"
        @test s1.spectrumType == :centroid
        @test s1.tic ≈ 4800.0   # sum of intensities
        @test s1.basePeakMz ≈ 150.0
        @test s1.basePeakIntensity ≈ 2000.0
        @test length(s1.mz) == 4
        @test s1.mz ≈ [110.0, 150.0, 200.0, 250.0]
        @test s1.int ≈ [500.0, 2000.0, 1500.0, 800.0]
        @test s1.metadata["title"] == "Spectrum 1, scan 101"

        # Scan 2
        s2 = scans[2]
        @test s2.num == 2
        @test s2.rt ≈ 1.0
        @test s2.precursor ≈ 500.0
        @test s2.chargeState == 3

        # Scan 3
        s3 = scans[3]
        @test s3.num == 3
        @test s3.precursor ≈ 600.0
        @test s3.chargeState == 1
        @test length(s3.mz) == 4

        # retention_time
        rt = MassJ.retention_time("test.mgf")
        @test length(rt) == 3
        @test rt ≈ [0.5, 1.0, 1.5]

        # chromatogram
        chrom = MassJ.chromatogram("test.mgf")
        @test length(chrom.rt) == 3
        @test chrom.maxic ≈ 4800.0

        # extract (filter by precursor)
        sub = MassJ.extract("test.mgf", MassJ.Precursor(500.0))
        @test length(sub) == 1
        @test sub[1].precursor ≈ 500.0

        # average
        avg = MassJ.average("test.mgf")
        @test avg isa MassJ.MSscans

        # New MSscan fields
        @test s1.mobilityType == :none
        @test s1.driftTime ≈ -1.0
        @test s1.compensationVoltage ≈ 0.0
    end
end


function test_msp()
    @testset "MSP format" begin
        # info
        inf = MassJ.info("test.msp")
        @test inf[1] == "3 scans"
        @test any(contains(s, "MS2+") for s in inf)
        @test any(contains(s, "MS2-") for s in inf)
        @test any(contains(s, "MS1+") for s in inf)

        inf_v = MassJ.info("test.msp", verbose=true)
        @test any(contains(s, "Caffeine") for s in inf_v)
        @test any(contains(s, "Aspirin") for s in inf_v)

        # load all
        scans = MassJ.load("test.msp")
        @test length(scans) == 3
        @test eltype(scans) == MassJ.MSscan

        # Scan 1: Caffeine MS2+
        s1 = scans[1]
        @test s1.num == 1
        @test s1.level == 2
        @test s1.polarity == "+"
        @test s1.precursor ≈ 195.0877
        @test s1.collisionEnergy ≈ 20.0
        @test s1.spectrumType == :centroid
        @test length(s1.mz) == 5
        @test s1.mz ≈ [42.0338, 69.0447, 110.0713, 138.0662, 195.0877]
        @test s1.int ≈ [1200.0, 3500.0, 8900.0, 45000.0, 120000.0]
        @test s1.basePeakMz ≈ 195.0877
        @test s1.basePeakIntensity ≈ 120000.0
        @test s1.tic ≈ sum(s1.int)
        @test s1.metadata["name"] == "Caffeine"
        @test s1.metadata["formula"] == "C8H10N4O2"
        @test s1.metadata["inchikey"] == "RYYVLZVUVIJVGH-UHFFFAOYSA-N"
        @test s1.metadata["db_id"] == "MSP001"
        @test s1.metadata["comments"] == "test spectrum 1"

        # Scan 2: Aspirin MS2-
        s2 = scans[2]
        @test s2.num == 2
        @test s2.level == 2
        @test s2.polarity == "-"
        @test s2.precursor ≈ 179.0344
        @test s2.collisionEnergy ≈ 15.0
        @test length(s2.mz) == 3
        @test s2.metadata["name"] == "Aspirin"
        @test s2.metadata["cas"] == "50-78-2"

        # Scan 3: Glucose MS1+ (semicolon format)
        s3 = scans[3]
        @test s3.num == 3
        @test s3.level == 1
        @test s3.polarity == "+"
        @test s3.precursor ≈ 0.0
        @test length(s3.mz) == 4
        @test s3.mz ≈ [180.0634, 181.0668, 182.0701, 183.0735]
        @test s3.int ≈ [8000.0, 900.0, 50.0, 5.0]
        @test s3.metadata["name"] == "Glucose"

        # retention_time
        rt = MassJ.retention_time("test.msp")
        @test length(rt) == 3
        @test all(rt .≈ 0.0)  # no RT in MSP test file

        # chromatogram
        chrom = MassJ.chromatogram("test.msp")
        @test length(chrom.rt) == 3

        # extract
        ms2_pos = MassJ.extract("test.msp", MassJ.Polarity("+"), MassJ.Level(2))
        @test length(ms2_pos) == 1
        @test ms2_pos[1].metadata["name"] == "Caffeine"

        ms2_neg = MassJ.extract("test.msp", MassJ.Polarity("-"))
        @test length(ms2_neg) == 1
        @test ms2_neg[1].metadata["name"] == "Aspirin"

        # average
        avg = MassJ.average("test.msp", MassJ.Level(2))
        @test avg isa MassJ.MSscans
    end
end


function test_imzml()
    @testset "imzML format" begin
        # info
        inf = MassJ.info("test.imzML")
        @test any(contains(s, "4 spectra") for s in inf)
        @test any(contains(s, "MS1+") for s in inf)

        inf_v = MassJ.info("test.imzML", verbose=true)
        @test any(contains(s, "processed") for s in inf_v)
        @test any(contains(s, "2 x 2") for s in inf_v)

        # load all
        scans = MassJ.load("test.imzML")
        @test length(scans) == 4
        @test eltype(scans) == MassJ.MSscan

        # Scan 1: position (1,1)
        s1 = scans[1]
        @test s1.num == 1
        @test s1.level == 1
        @test s1.polarity == "+"
        @test s1.spectrumType == :profile
        @test s1.tic ≈ 8000.0
        @test length(s1.mz) == 3
        @test s1.mz ≈ [100.0, 200.0, 300.0]
        @test s1.int ≈ [1000.0, 5000.0, 2000.0]
        @test s1.metadata["position_x"] == 1
        @test s1.metadata["position_y"] == 1

        # Scan 2: position (2,1)
        s2 = scans[2]
        @test s2.num == 2
        @test s2.tic ≈ 7800.0
        @test s2.mz ≈ [100.0, 200.0, 300.0]
        @test s2.int ≈ [800.0, 3000.0, 4000.0]
        @test s2.metadata["position_x"] == 2
        @test s2.metadata["position_y"] == 1

        # Scan 3: position (1,2), 2 peaks
        s3 = scans[3]
        @test s3.num == 3
        @test length(s3.mz) == 2
        @test s3.mz ≈ [150.0, 250.0]
        @test s3.int ≈ [2000.0, 6000.0]
        @test s3.metadata["position_x"] == 1
        @test s3.metadata["position_y"] == 2

        # Scan 4: position (2,2)
        s4 = scans[4]
        @test s4.num == 4
        @test length(s4.mz) == 3
        @test s4.mz ≈ [150.0, 250.0, 350.0]
        @test s4.int ≈ [1500.0, 4000.0, 3000.0]
        @test s4.metadata["position_x"] == 2
        @test s4.metadata["position_y"] == 2

        # retention_time
        rt = MassJ.retention_time("test.imzML")
        @test length(rt) == 4

        # chromatogram
        chrom = MassJ.chromatogram("test.imzML")
        @test length(chrom.rt) == 4

        # extract
        all_scans = MassJ.extract("test.imzML", MassJ.Level(1))
        @test length(all_scans) == 4

        # average
        avg = MassJ.average("test.imzML")
        @test avg isa MassJ.MSscans
    end
end


function test_composed_predicates()
    scans = MassJ.load("test.mzXML")

    @testset "Composed predicates - empty filter returns all scans" begin
        sub = MassJ.extract(scans)
        @test length(sub) == length(scans)

        chrom = MassJ.chromatogram(scans)
        @test length(chrom.rt) == length(scans)
    end

    @testset "Composed predicates - no match returns ErrorException" begin
        @test MassJ.extract(scans, MassJ.Level(99))       isa ErrorException
        @test MassJ.chromatogram(scans, MassJ.Level(99))  isa ErrorException
        @test MassJ.average(scans, MassJ.Level(99))       isa ErrorException
    end

    @testset "Composed predicates - single match in average returns MSscan" begin
        result = MassJ.average(scans, MassJ.Scan(1))
        @test result isa MassJ.MSscan
        @test result.num == 1
    end

    @testset "Composed predicates - AND semantics equivalence" begin
        # Single-pass composition must yield the same scans as stepwise filtering.
        combined = MassJ.extract(scans, MassJ.Level(2), MassJ.Polarity("+"))
        stepwise = MassJ.extract(MassJ.extract(scans, MassJ.Level(2)), MassJ.Polarity("+"))
        @test [s.num for s in combined] == [s.num for s in stepwise]
    end

    @testset "Composed predicates - multiple disjoint RT ranges" begin
        # Two RT intervals that together cover scans at RT≈0.14, 0.73 and 4.34 in test.mzXML.
        ms = MassJ.average(scans, MassJ.RT([[0.0, 1.0], [4.0, 5.0]]), stats = false)
        @test ms isa MassJ.MSscans
    end
end


function test_yields()
    @testset "Yields - integrate_window, yields, normalize_*" begin

        # integrate_window on synthetic arrays: triangle from 0..1..0 over m/z 0..2
        mz   = [0.0, 1.0, 2.0]
        int  = [0.0, 1.0, 0.0]
        @test MassJ.integrate_window(mz, int, 0.0, 2.0) ≈ 1.0   # 2 * 1 / 2
        @test MassJ.integrate_window(mz, int, 0.0, 1.0) ≈ 0.5
        @test MassJ.integrate_window(mz, int, 5.0, 6.0) == 0.0  # no points
        @test MassJ.integrate_window(mz, int, 1.0, 1.0) == 0.0  # 1 point only
        @test MassJ.integrate_window(mz, int, 2.0, 0.0) ≈ 1.0   # swapped bounds

        # integrate_window on a real MSscan
        scans = MassJ.load("test.mzXML")
        a = MassJ.integrate_window(scans[1], 400.0, 500.0)
        @test a > 0 && isfinite(a)

        # Peak constructor swaps when mz1 > mz2
        p = MassJ.Peak(200.0, 100.0, "swap")
        @test p.mz1 == 100.0 && p.mz2 == 200.0 && p.label == "swap"

        # yields(files, peaks; x) — two-file series using the same fixture
        peaks = [MassJ.Peak(400.0, 500.0, "low"),
                 MassJ.Peak(800.0, 900.0, "high")]
        yc = MassJ.yields(["test.mzXML", "test.mzXML"], peaks;
                          x = [3.5, 4.0], xlabel = "photon energy (eV)")
        @test yc isa MassJ.YieldCurve
        @test size(yc.yields)    == (2, 2)
        @test size(yc.found_mz)  == (2, 2)
        @test all(isnan, yc.found_mz)         # all fixed Peak → all NaN
        @test yc.x               == [3.5, 4.0]
        @test yc.xlabel          == "photon energy (eV)"
        @test yc.labels          == ["low", "high"]
        @test yc.windows         == [(400.0, 500.0), (800.0, 900.0)]
        @test yc.tic[1]          ≈ yc.yields[1, 1] + yc.yields[1, 2]
        @test yc.yields[1, :]    ≈ yc.yields[2, :]    # same source file twice

        # length mismatch
        @test_throws ErrorException MassJ.yields(["test.mzXML"], peaks; x = [1.0, 2.0])

        # read_peaklist — round-trip through a temp CSV (with header)
        tmp = tempname() * ".csv"
        open(tmp, "w") do io
            write(io, "mz1,mz2,label\n")
            write(io, "400.0,500.0,low\n")
            write(io, "800.0,900.0,high\n")
        end
        pl = MassJ.read_peaklist(tmp)
        @test length(pl) == 2
        @test pl[1].mz1 == 400.0 && pl[1].mz2 == 500.0 && pl[1].label == "low"
        rm(tmp)

        # read_peaklist — no header
        tmp2 = tempname() * ".csv"
        open(tmp2, "w") do io
            write(io, "100.0,200.0,A\n")
        end
        pl2 = MassJ.read_peaklist(tmp2)
        @test length(pl2) == 1 && pl2[1].label == "A"
        rm(tmp2)

        # normalize_tic: rows of peak columns sum to 1
        yn = MassJ.normalize_tic(yc)
        @test yn isa MassJ.YieldCurve
        @test sum(yn.yields[1, :]) ≈ 1.0
        @test sum(yn.yields[2, :]) ≈ 1.0
        @test yn.tic == yc.tic   # raw totals preserved
        @test yn.metadata["normalize_tic"] == true

        # normalize_flux: divide by a constant flux of 2.0
        # — `#`-prefixed header lines are stripped as comments
        fluxpath = tempname() * ".txt"
        open(fluxpath, "w") do io
            write(io, "# header line 1\n")
            write(io, "# header line 2\n")
            write(io, "3.0  2.0\n")
            write(io, "5.0  2.0\n")
        end
        yf = MassJ.normalize_flux(yc, fluxpath)
        @test yf.yields ≈ yc.yields ./ 2.0
        @test yf.tic    ≈ yc.tic    ./ 2.0
        @test yf.metadata["normalize_flux"] == fluxpath
        rm(fluxpath)

        # normalize_flux: text header (no #) auto-detected and skipped
        flux_text = tempname() * ".txt"
        open(flux_text, "w") do io
            write(io, "energy flux\n")     # text header
            write(io, "==== ====\n")       # decorative
            write(io, "3.0  2.0\n")
            write(io, "5.0  2.0\n")
        end
        yf_t = MassJ.normalize_flux(yc, flux_text)
        @test yf_t.yields ≈ yc.yields ./ 2.0
        rm(flux_text)

        # normalize_flux: # comments mixed with data anywhere in the file
        flux_mixed = tempname() * ".txt"
        open(flux_mixed, "w") do io
            write(io, "# preamble\n")
            write(io, "energy   flux\n")   # text header below comments
            write(io, "3.0  2.0\n")
            write(io, "# mid-file note\n")
            write(io, "5.0  2.0  # trailing comment\n")
        end
        yf_m = MassJ.normalize_flux(yc, flux_mixed)
        @test yf_m.yields ≈ yc.yields ./ 2.0
        rm(flux_mixed)

        # normalize_flux: explicit skipstart override
        # — top row is "1 1" (looks numeric, but is a unit/scale row to skip)
        flux_skip = tempname() * ".txt"
        open(flux_skip, "w") do io
            write(io, "1 1\n")             # would otherwise be parsed as data
            write(io, "3.0  2.0\n")
            write(io, "5.0  2.0\n")
        end
        yf_s = MassJ.normalize_flux(yc, flux_skip; skipstart = 1)
        @test yf_s.yields ≈ yc.yields ./ 2.0
        rm(flux_skip)

        # write_csv round-trip — header + correct row count
        outpath = tempname() * ".csv"
        MassJ.write_csv(yc, outpath)
        lines = readlines(outpath)
        @test length(lines) == 3                                   # header + 2 rows
        @test lines[1] == "photon energy (eV),low,high,TIC"
        rm(outpath)

        # plot recipe smoke test
        @test typeof(plot(yc)) == Plots.Plot{Plots.GRBackend}
    end
end


function test_yields_targetpeak()
    @testset "TargetPeak resolution methods + Peak(mz, label; tol/ppm)" begin

        # Peak(mz, label; tol) — eager fixed window from a single m/z
        p = MassJ.Peak(100.0, "x"; tol = 0.5)
        @test p.mz1 == 99.5 && p.mz2 == 100.5 && p.label == "x"

        # Peak(mz, label; ppm)
        p2 = MassJ.Peak(1000.0, "y"; ppm = 5.0)
        @test p2.mz1 ≈ 1000.0 - 1000.0 * 5e-6
        @test p2.mz2 ≈ 1000.0 + 1000.0 * 5e-6

        # tol/ppm mutual exclusion
        @test_throws ErrorException MassJ.Peak(100.0, "z"; tol = 0.5, ppm = 5.0)
        @test_throws ErrorException MassJ.Peak(100.0, "z")

        # TargetPeak defaults and method validation
        tp = MassJ.TargetPeak(100.0, "a"; tol = 0.5)
        @test tp.mz == 100.0 && tp.tol == 0.5
        @test tp.method === :local_max && tp.edges == 0.1

        tp2 = MassJ.TargetPeak(100.0, "b"; tol = 0.5, method = :edges, edges = 0.2)
        @test tp2.method === :edges && tp2.edges == 0.2

        @test_throws ErrorException MassJ.TargetPeak(100.0, "c"; tol = 0.5, method = :bad)
        @test_throws ErrorException MassJ.TargetPeak(100.0, "c")     # neither tol nor ppm

        # Locate the global max in the averaged test fixture
        spec = MassJ.average("test.mzXML")
        peak_idx  = argmax(spec.int)
        target_mz = spec.mz[peak_idx]

        # :local_max — should snap exactly onto the sample-grid maximum
        peaks_lm = [MassJ.TargetPeak(target_mz - 0.05, "lm"; tol = 0.2)]
        yc_lm    = MassJ.yields(["test.mzXML"], peaks_lm; x = [1.0])
        @test yc_lm.found_mz[1, 1] ≈ target_mz
        @test yc_lm.yields[1, 1]   > 0

        # :edges — same location, possibly different window width
        peaks_ed = [MassJ.TargetPeak(target_mz - 0.05, "ed";
                                      tol = 0.2, method = :edges)]
        yc_ed    = MassJ.yields(["test.mzXML"], peaks_ed; x = [1.0])
        @test yc_ed.found_mz[1, 1] ≈ target_mz
        @test yc_ed.yields[1, 1]   > 0

        # :centroid — uses the package's centroid(). Wider tol + TBPD so the
        # search window definitely contains at least one centroid in this fixture.
        peaks_cn = [MassJ.TargetPeak(target_mz, "cn";
                                      tol = 5.0, method = :centroid)]
        yc_cn    = MassJ.yields(["test.mzXML"], peaks_cn; x = [1.0],
                                centroid_method = MassJ.TBPD(:gauss, 4500., 0.2))
        @test isfinite(yc_cn.found_mz[1, 1])
        @test abs(yc_cn.found_mz[1, 1] - target_mz) < 5.0
        @test yc_cn.yields[1, 1] > 0

        # Mixed peak list — Peak gives NaN, TargetPeak gives a located m/z
        peaks_mx = [MassJ.Peak(400.0, 500.0, "static"),
                    MassJ.TargetPeak(target_mz, "lazy"; tol = 0.2)]
        yc_mx    = MassJ.yields(["test.mzXML"], peaks_mx; x = [1.0])
        @test isnan(yc_mx.found_mz[1, 1])
        @test yc_mx.found_mz[1, 2] ≈ target_mz

        # read_peaklist 2-col → TargetPeak with kwarg defaults
        tmp2 = tempname() * ".csv"
        open(tmp2, "w") do io
            write(io, "mz,label\n")
            write(io, "100.0,A\n")
            write(io, "200.0,B\n")
        end
        pl2 = MassJ.read_peaklist(tmp2; tol = 0.3, method = :edges)
        @test length(pl2) == 2
        @test pl2[1] isa MassJ.TargetPeak
        @test pl2[1].mz == 100.0 && pl2[1].tol == 0.3 && pl2[1].method === :edges
        rm(tmp2)

        # read_peaklist 4-col → TargetPeak with per-row tol + method
        tmp4 = tempname() * ".csv"
        open(tmp4, "w") do io
            write(io, "mz,tol,method,label\n")
            write(io, "100.0,0.3,local_max,A\n")
            write(io, "200.0,0.5,edges,B\n")
        end
        pl4 = MassJ.read_peaklist(tmp4)
        @test length(pl4) == 2
        @test pl4[1].method === :local_max && pl4[1].tol == 0.3
        @test pl4[2].method === :edges     && pl4[2].tol == 0.5
        rm(tmp4)

        # normalize_tic preserves found_mz
        yn = MassJ.normalize_tic(yc_lm)
        @test yn.found_mz == yc_lm.found_mz

        # drop_peaks — single string
        peaks_two = [MassJ.Peak(400.0, 500.0, "low"),
                     MassJ.Peak(800.0, 900.0, "high")]
        yc2 = MassJ.yields(["test.mzXML", "test.mzXML"], peaks_two;
                           x = [1.0, 2.0])
        d1 = MassJ.drop_peaks(yc2, "low")
        @test d1.labels        == ["high"]
        @test size(d1.yields)  == (2, 1)
        @test d1.yields[:, 1]  ≈  yc2.yields[:, 2]
        @test d1.windows       == [(800.0, 900.0)]
        @test size(d1.found_mz) == (2, 1)
        @test d1.tic           == yc2.tic         # tic unchanged by design

        # drop_peaks — vector of labels
        d2 = MassJ.drop_peaks(yc2, ["high"])
        @test d2.labels == ["low"]

        # drop_peaks — drop everything
        d3 = MassJ.drop_peaks(yc2, ["low", "high"])
        @test isempty(d3.labels)
        @test size(d3.yields) == (2, 0)

        # drop_peaks — missing labels silently ignored
        d4 = MassJ.drop_peaks(yc2, ["nonexistent"])
        @test d4.labels == yc2.labels
        @test d4.yields == yc2.yields
    end
end


function test_yields_errors()
    @testset "YieldCurve error propagation" begin

        # MSscans from `average("test.mzXML")` carries variance over 6 scans,
        # so yields_err and tic_err should all be finite.
        peaks = [MassJ.Peak(400.0, 500.0, "low"),
                 MassJ.Peak(800.0, 900.0, "high")]
        yc = MassJ.yields(["test.mzXML", "test.mzXML"], peaks; x = [1.0, 2.0])

        @test size(yc.yields_err) == (2, 2)
        @test all(isfinite, yc.yields_err)
        @test all(yc.yields_err .>= 0)
        @test size(yc.tic_err) == (2,)
        @test all(isfinite, yc.tic_err)
        for i in 1:2
            @test yc.tic_err[i] ≈ sqrt(sum(abs2, yc.yields_err[i, :]))
        end

        # MSscan path: no variance available → NaN error
        scans = MassJ.load("test.mzXML")
        _, σ_one = MassJ._integrate_window_with_err(scans[1], 400.0, 500.0)
        @test isnan(σ_one)

        # normalize_tic: standard division error propagation
        yn = MassJ.normalize_tic(yc)
        @test all(isfinite, yn.yields_err)
        for i in 1:2, p in 1:2
            σ_y = yc.yields_err[i, p]
            y   = yc.yields[i, p]
            t   = yc.tic[i]
            σ_t = yc.tic_err[i]
            expected = sqrt((σ_y / t)^2 + (y * σ_t / (t * t))^2)
            @test yn.yields_err[i, p] ≈ expected
        end
        @test yn.tic     == yc.tic           # raw totals preserved
        @test yn.tic_err == yc.tic_err

        # normalize_flux: 2-col file → 10% default σ_φ
        fluxpath = tempname() * ".txt"
        open(fluxpath, "w") do io
            write(io, "# flux 10%\n")
            write(io, "1.0  2.0\n")
            write(io, "3.0  2.0\n")
        end
        yf = MassJ.normalize_flux(yc, fluxpath)
        @test yf.metadata["normalize_flux_err_pct"] == 0.10
        for i in 1:2, p in 1:2
            φ, σφ = 2.0, 0.10 * 2.0
            σ_y   = yc.yields_err[i, p]
            y     = yc.yields[i, p]
            expected = sqrt((σ_y / φ)^2 + (y * σφ / (φ * φ))^2)
            @test yf.yields_err[i, p] ≈ expected
        end
        rm(fluxpath)

        # normalize_flux: custom flux_err_pct kwarg
        fluxpath2 = tempname() * ".txt"
        open(fluxpath2, "w") do io
            write(io, "1.0  2.0\n")
            write(io, "3.0  2.0\n")
        end
        yf5 = MassJ.normalize_flux(yc, fluxpath2; flux_err_pct = 0.05)
        @test yf5.metadata["normalize_flux_err_pct"] == 0.05
        rm(fluxpath2)

        # normalize_flux: 3-col file → σ_φ from the file
        flux3 = tempname() * ".txt"
        open(flux3, "w") do io
            write(io, "# x flux sigma\n")
            write(io, "1.0  2.0  0.1\n")
            write(io, "3.0  2.0  0.1\n")
        end
        yf3 = MassJ.normalize_flux(yc, flux3)
        for i in 1:2, p in 1:2
            φ, σφ = 2.0, 0.1
            σ_y   = yc.yields_err[i, p]
            y     = yc.yields[i, p]
            expected = sqrt((σ_y / φ)^2 + (y * σφ / (φ * φ))^2)
            @test yf3.yields_err[i, p] ≈ expected
        end
        rm(flux3)

        # normalize_flux: jagged / empty 3rd column should not throw
        # (readdlm sometimes pads files with an empty 3rd col when rows are
        #  uneven or have trailing whitespace) — per-row fallback to pct.
        # Use xf == yc.x so there's no interpolation between σ values.
        flux_jag = tempname() * ".txt"
        open(flux_jag, "w") do io
            write(io, "# x flux [sigma]\n")
            write(io, "1.0  2.0  0.1\n")   # σ present
            write(io, "2.0  2.0\n")        # σ missing → falls back to pct
        end
        yfj = MassJ.normalize_flux(yc, flux_jag; flux_err_pct = 0.10)
        @test all(isfinite, yfj.yields_err)
        # Row 1 (yc.x=1.0) lands on xf[1] → σφ = 0.1
        # Row 2 (yc.x=2.0) lands on xf[2] → σφ = 0.10·|2.0| = 0.2 (pct fallback)
        for p in 1:2
            φ = 2.0
            σ_y1 = yc.yields_err[1, p]; y1 = yc.yields[1, p]
            σ_y2 = yc.yields_err[2, p]; y2 = yc.yields[2, p]
            @test yfj.yields_err[1, p] ≈ sqrt((σ_y1 / φ)^2 + (y1 * 0.1 / (φ * φ))^2)
            @test yfj.yields_err[2, p] ≈ sqrt((σ_y2 / φ)^2 + (y2 * 0.2 / (φ * φ))^2)
        end
        rm(flux_jag)

        # normalize_flux: non-numeric 3rd column (e.g. date string in DESIRS
        # beamline log files) — each row's σ falls back to pct.
        flux_str = tempname() * ".txt"
        open(flux_str, "w") do io
            write(io, "# DESIRS-style header\n")
            write(io, "Energy (eV)  flux  timestamp\n")          # text header
            write(io, "1.0  2.0  Fri Apr 17 12:27:54 2026\n")
            write(io, "2.0  2.0  Fri Apr 17 12:28:02 2026\n")
        end
        yfs = MassJ.normalize_flux(yc, flux_str; flux_err_pct = 0.10)
        @test all(isfinite, yfs.yields_err)
        rm(flux_str)

        # drop_peaks slices yields_err; tic_err unchanged by design
        d = MassJ.drop_peaks(yc, "low")
        @test size(d.yields_err) == (2, 1)
        @test d.yields_err[:, 1] == yc.yields_err[:, 2]
        @test d.tic_err == yc.tic_err

        # Plot with ribbon still works
        @test typeof(plot(yc)) == Plots.Plot{Plots.GRBackend}
    end
end


tests()
test_isotopes()
test_deconvolution()
test_interpolation_import()
test_mzml()
test_mgf()
test_msp()
test_imzml()
test_composed_predicates()
test_yields()
test_yields_targetpeak()
test_yields_errors()
