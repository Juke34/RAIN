#!/usr/bin/env python3
"""
Tests pour les colonnes TotalSites, ObservedSites et QualifiedSites

Pour exécuter: 
    cd /workspaces/rain/bin
    PYTHONPATH=/workspaces/rain/bin python3 -m unittest pluviometer.test_site_counts -v
"""

import unittest
import numpy as np
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation

# Import du package pluviometer
from pluviometer.multi_counter import MultiCounter
from pluviometer.site_filter import SiteFilter
from pluviometer.utils import RNASiteVariantData
from pluviometer.__main__ import AggregatePositions


class TestAggregatePositions(unittest.TestCase):
    """Tests pour la classe AggregatePositions"""

    def test_initial_state(self):
        """Test l'état initial"""
        pos = AggregatePositions(strand=1)
        self.assertEqual(pos.start, float('inf'))
        self.assertEqual(pos.end, 0)
        self.assertEqual(pos.strand, 1)

    def test_update_from_simple_feature(self):
        """Test mise à jour avec une SimpleLocation"""
        pos = AggregatePositions(strand=1)
        feature = SeqFeature(location=SimpleLocation(10, 20), id="test", type="exon")
        pos.update_from_feature(feature)
        
        self.assertEqual(pos.start, 10)
        self.assertEqual(pos.end, 20)
        self.assertEqual(pos.strand, 1)

    def test_update_from_multiple_features(self):
        """Test mise à jour avec plusieurs features"""
        pos = AggregatePositions(strand=1)
        
        feature1 = SeqFeature(location=SimpleLocation(10, 20), id="exon1", type="exon")
        feature2 = SeqFeature(location=SimpleLocation(5, 15), id="exon2", type="exon")
        feature3 = SeqFeature(location=SimpleLocation(30, 40), id="exon3", type="exon")
        
        pos.update_from_feature(feature1)
        pos.update_from_feature(feature2)
        pos.update_from_feature(feature3)
        
        # Start doit être le minimum (5), end le maximum (40)
        self.assertEqual(pos.start, 5)
        self.assertEqual(pos.end, 40)

    def test_merge_positions(self):
        """Test fusion de deux AggregatePositions"""
        pos1 = AggregatePositions(strand=1)
        pos1.start = 10
        pos1.end = 20
        
        pos2 = AggregatePositions(strand=1)
        pos2.start = 5
        pos2.end = 30
        
        pos1.merge(pos2)
        
        self.assertEqual(pos1.start, 5)
        self.assertEqual(pos1.end, 30)

    def test_to_strings_valid(self):
        """Test conversion en strings avec positions valides"""
        pos = AggregatePositions(strand=1)
        pos.start = 10
        pos.end = 20
        
        start_str, end_str, strand_str = pos.to_strings()
        
        # Start converti en 1-based
        self.assertEqual(start_str, "11")
        self.assertEqual(end_str, "20")
        self.assertEqual(strand_str, "1")

    def test_to_strings_empty(self):
        """Test conversion en strings avec positions non initialisées"""
        pos = AggregatePositions(strand=0)
        
        start_str, end_str, strand_str = pos.to_strings()
        
        self.assertEqual(start_str, ".")
        self.assertEqual(end_str, ".")
        self.assertEqual(strand_str, ".")


class TestMultiCounterSiteCounts(unittest.TestCase):
    """Tests pour les compteurs de sites dans MultiCounter"""

    def test_site_base_pairings_counts_sites_not_reads(self):
        """Test que SiteBasePairingsObserved compte les SITES (pas les reads) avec chaque pairing"""
        # Reproduit le cas du test minimal:
        # pos1: A ref, cov=10, [10,0,0,0] -> AA pairing passes (10>=5)
        # pos2: C ref, cov=12, [0,9,0,3]  -> CC passes (9>=5), CT filtered (3<5)
        # pos4: T ref, cov=2              -> tout filtré (cov<5)
        filter = SiteFilter(cov_threshold=5, edit_threshold=5)
        counter = MultiCounter(filter)

        counter.update(RNASiteVariantData("21", 0, 0, 1, 10, 37.0, np.array([10,0,0,0,0], dtype=np.int32), 0.0))  # A site
        counter.update(RNASiteVariantData("21", 1, 1, 1, 12, 37.0, np.array([0,9,0,3,0], dtype=np.int32), 0.0))  # C site
        counter.update(RNASiteVariantData("21", 3, 3, 1,  2, 37.0, np.array([0,0,0,2,0], dtype=np.int32), 0.0))  # T site, cov<5

        # SiteBasePairingsObserved: 1 site AA, 1 site CC, 0 CT (filtré), 0 TT (cov<5)
        self.assertEqual(counter.edit_site_freqs[0, 0], 1)   # AA: 1 site
        self.assertEqual(counter.edit_site_freqs[0, 1], 0)   # AC: 0
        self.assertEqual(counter.edit_site_freqs[1, 1], 1)   # CC: 1 site
        self.assertEqual(counter.edit_site_freqs[1, 3], 0)   # CT: 0 (3 reads < edit_threshold)
        self.assertEqual(counter.edit_site_freqs[3, 3], 0)   # TT: 0 (cov < cov_threshold)

        # ReadBasePairings: compte les reads bruts (non filtrés)
        self.assertEqual(counter.edit_read_freqs[0, 0], 10)  # AA: 10 reads
        self.assertEqual(counter.edit_read_freqs[1, 1], 9)   # CC: 9 reads
        self.assertEqual(counter.edit_read_freqs[1, 3], 3)   # CT: 3 reads (non filtrés dans ReadBasePairings)
        self.assertEqual(counter.edit_read_freqs[3, 3], 2)   # TT: 2 reads

    def test_observed_sites_count(self):
        """Test que ObservedSites compte toutes les observations"""
        filter = SiteFilter(cov_threshold=5, edit_threshold=2)
        counter = MultiCounter(filter)
        
        # Ajouter 3 sites avec différentes couvertures
        site1 = RNASiteVariantData(
            seqid="chr1", position=10, reference=0,  # A
            strand=1, coverage=10, mean_quality=30.0,
            frequencies=np.array([8, 2, 0, 0, 0], dtype=np.int32),
            score=0.0
        )
        site2 = RNASiteVariantData(
            seqid="chr1", position=20, reference=1,  # C
            strand=1, coverage=3, mean_quality=30.0,
            frequencies=np.array([1, 2, 0, 0, 0], dtype=np.int32),
            score=0.0
        )
        site3 = RNASiteVariantData(
            seqid="chr1", position=30, reference=2,  # G
            strand=1, coverage=8, mean_quality=30.0,
            frequencies=np.array([0, 3, 5, 0, 0], dtype=np.int32),
            score=0.0
        )
        
        counter.update(site1)
        counter.update(site2)
        counter.update(site3)
        
        # ObservedSites = sum de genome_base_freqs = 3 sites
        observed = counter.genome_base_freqs.sum()
        self.assertEqual(observed, 3)
        
        # Vérifier la distribution par base
        self.assertEqual(counter.genome_base_freqs[0], 1)  # A
        self.assertEqual(counter.genome_base_freqs[1], 1)  # C
        self.assertEqual(counter.genome_base_freqs[2], 1)  # G

    def test_qualified_sites_count(self):
        """Test que QualifiedSites compte seulement les sites passant le filtre de couverture"""
        filter = SiteFilter(cov_threshold=5, edit_threshold=2)
        counter = MultiCounter(filter)
        
        # Site 1: couverture 10 >= 5 → qualifié (référence A)
        site1 = RNASiteVariantData(
            seqid="chr1", position=10, reference=0,  # A
            strand=1, coverage=10, mean_quality=30.0,
            frequencies=np.array([8, 2, 0, 0, 0], dtype=np.int32),
            score=0.0
        )
        # Site 2: couverture 3 < 5 → non qualifié (référence C)
        site2 = RNASiteVariantData(
            seqid="chr1", position=20, reference=1,  # C
            strand=1, coverage=3, mean_quality=30.0,
            frequencies=np.array([1, 2, 0, 0, 0], dtype=np.int32),
            score=0.0
        )
        # Site 3: couverture 8 >= 5 → qualifié (référence G)
        site3 = RNASiteVariantData(
            seqid="chr1", position=30, reference=2,  # G
            strand=1, coverage=8, mean_quality=30.0,
            frequencies=np.array([0, 3, 5, 0, 0], dtype=np.int32),
            score=0.0
        )
        
        counter.update(site1)
        counter.update(site2)
        counter.update(site3)
        
        # QualifiedSites = 2 (site1 et site3)
        self.assertEqual(counter.filtered_sites_count, 2)
        
        # ObservedSites = 3 (tous)
        self.assertEqual(counter.genome_base_freqs.sum(), 3)
        
        # QualifiedBases: seulement A et G (site1 et site3)
        self.assertEqual(counter.filtered_base_freqs[0], 1)  # A
        self.assertEqual(counter.filtered_base_freqs[1], 0)  # C (non qualifié)
        self.assertEqual(counter.filtered_base_freqs[2], 1)  # G
        self.assertEqual(counter.filtered_base_freqs[3], 0)  # T

    def test_merge_preserves_counts(self):
        """Test que merge préserve les compteurs de sites et de bases"""
        filter = SiteFilter(cov_threshold=5, edit_threshold=2)
        
        counter1 = MultiCounter(filter)
        site1 = RNASiteVariantData(
            seqid="chr1", position=10, reference=0,  # A, coverage 10 >= 5
            strand=1, coverage=10, mean_quality=30.0,
            frequencies=np.array([8, 2, 0, 0, 0], dtype=np.int32),
            score=0.0
        )
        counter1.update(site1)
        
        counter2 = MultiCounter(filter)
        site2 = RNASiteVariantData(
            seqid="chr1", position=20, reference=1,  # C, coverage 3 < 5
            strand=1, coverage=3, mean_quality=30.0,
            frequencies=np.array([1, 2, 0, 0, 0], dtype=np.int32),
            score=0.0
        )
        site3 = RNASiteVariantData(
            seqid="chr1", position=30, reference=0,  # A, coverage 6 >= 5
            strand=1, coverage=6, mean_quality=30.0,
            frequencies=np.array([5, 1, 0, 0, 0], dtype=np.int32),
            score=0.0
        )
        counter2.update(site2)
        counter2.update(site3)
        
        # Merge counter2 into counter1
        counter1.merge(counter2)
        
        # Vérifier les totaux
        self.assertEqual(counter1.genome_base_freqs.sum(), 3)  # 1 + 2 = 3 sites observés
        self.assertEqual(counter1.filtered_sites_count, 2)  # 1 + 1 = 2 sites qualifiés (site1 et site3)
        
        # Vérifier les bases observées
        self.assertEqual(counter1.genome_base_freqs[0], 2)  # A: site1 + site3
        self.assertEqual(counter1.genome_base_freqs[1], 1)  # C: site2
        
        # Vérifier les bases qualifiées
        self.assertEqual(counter1.filtered_base_freqs[0], 2)  # A: site1 + site3 (tous deux qualifiés)
        self.assertEqual(counter1.filtered_base_freqs[1], 0)  # C: site2 non qualifié
        self.assertEqual(counter1.filtered_base_freqs.sum(), 2)  # Total = 2

    def test_qualified_bases_distribution(self):
        """Test que QualifiedBases reflète correctement la distribution des bases qualifiées"""
        filter = SiteFilter(cov_threshold=5, edit_threshold=2)
        counter = MultiCounter(filter)
        
        # Ajouter plusieurs sites avec différentes bases et couvertures
        sites = [
            (0, 10),  # A, qualifié
            (0, 6),   # A, qualifié
            (1, 3),   # C, non qualifié
            (1, 7),   # C, qualifié
            (2, 2),   # G, non qualifié
            (2, 5),   # G, qualifié
            (3, 8),   # T, qualifié
        ]
        
        for ref, cov in sites:
            site = RNASiteVariantData(
                seqid="chr1",
                position=0,
                reference=ref,
                strand=1,
                coverage=cov,
                mean_quality=30.0,
                frequencies=np.array([1, 0, 0, 0, 0], dtype=np.int32),
                score=0.0
            )
            counter.update(site)
        
        # ObservedBases: toutes les occurrences
        self.assertEqual(counter.genome_base_freqs[0], 2)  # A: 2
        self.assertEqual(counter.genome_base_freqs[1], 2)  # C: 2
        self.assertEqual(counter.genome_base_freqs[2], 2)  # G: 2
        self.assertEqual(counter.genome_base_freqs[3], 1)  # T: 1
        self.assertEqual(counter.genome_base_freqs.sum(), 7)
        
        # QualifiedBases: seulement avec cov >= 5
        self.assertEqual(counter.filtered_base_freqs[0], 2)  # A: 2 qualifiés (10, 6)
        self.assertEqual(counter.filtered_base_freqs[1], 1)  # C: 1 qualifié (7)
        self.assertEqual(counter.filtered_base_freqs[2], 1)  # G: 1 qualifié (5)
        self.assertEqual(counter.filtered_base_freqs[3], 1)  # T: 1 qualifié (8)
        self.assertEqual(counter.filtered_base_freqs.sum(), 5)
        self.assertEqual(counter.filtered_sites_count, 5)


class TestFeatureTotalSites(unittest.TestCase):
    """Tests pour le calcul de TotalSites des features"""

    def test_simple_feature_total_sites(self):
        """Test TotalSites pour une feature simple"""
        feature = SeqFeature(
            location=SimpleLocation(0, 100),
            id="exon1",
            type="exon"
        )
        
        # TotalSites = len(location)
        total_sites = len(feature.location)
        self.assertEqual(total_sites, 100)

    def test_compound_feature_total_sites(self):
        """Test TotalSites pour une feature avec CompoundLocation (chimaera)"""
        # Deux exons: [0-50] et [100-150]
        location = CompoundLocation([
            SimpleLocation(0, 50),
            SimpleLocation(100, 150)
        ])
        feature = SeqFeature(
            location=location,
            id="chimaera",
            type="exon"
        )
        
        # TotalSites = somme des longueurs sans gaps
        total_sites = len(feature.location)
        self.assertEqual(total_sites, 100)  # 50 + 50

    def test_gene_total_sites(self):
        """Test TotalSites pour un gene (span complet incluant introns)"""
        # Gene de 1 à 1000
        feature = SeqFeature(
            location=SimpleLocation(0, 1000),
            id="gene1",
            type="gene"
        )
        
        # TotalSites = span complet
        total_sites = len(feature.location)
        self.assertEqual(total_sites, 1000)


class TestAggregateTotalSites(unittest.TestCase):
    """Tests pour le calcul de TotalSites des aggregates"""

    def test_aggregate_from_continuous_features(self):
        """Test TotalSites pour aggregate de features contiguës"""
        pos = AggregatePositions(strand=1)
        
        # Trois exons adjacents
        exon1 = SeqFeature(location=SimpleLocation(0, 100), id="exon1", type="exon")
        exon2 = SeqFeature(location=SimpleLocation(100, 200), id="exon2", type="exon")
        exon3 = SeqFeature(location=SimpleLocation(200, 300), id="exon3", type="exon")
        
        pos.update_from_feature(exon1)
        pos.update_from_feature(exon2)
        pos.update_from_feature(exon3)
        
        # TotalSites = end - start = 300 - 0 = 300
        total_sites = pos.end - pos.start
        self.assertEqual(total_sites, 300)

    def test_aggregate_from_separated_features(self):
        """Test TotalSites pour aggregate de features séparées"""
        pos = AggregatePositions(strand=1)
        
        # Deux exons avec un gap
        exon1 = SeqFeature(location=SimpleLocation(0, 100), id="exon1", type="exon")
        exon2 = SeqFeature(location=SimpleLocation(200, 300), id="exon2", type="exon")
        
        pos.update_from_feature(exon1)
        pos.update_from_feature(exon2)
        
        # TotalSites = end - start = 300 - 0 = 300 (inclut le gap)
        total_sites = pos.end - pos.start
        self.assertEqual(total_sites, 300)

    def test_aggregate_empty(self):
        """Test TotalSites pour aggregate vide"""
        pos = AggregatePositions(strand=0)
        
        # Aucune feature ajoutée
        start_str, end_str, strand_str = pos.to_strings()
        
        # Devrait retourner "."
        self.assertEqual(start_str, ".")
        self.assertEqual(end_str, ".")


class TestIntegrationSiteCounts(unittest.TestCase):
    """Tests d'intégration pour les trois colonnes"""

    def test_all_three_columns(self):
        """Test que les trois colonnes sont cohérentes"""
        filter = SiteFilter(cov_threshold=5, edit_threshold=2)
        counter = MultiCounter(filter)
        
        # 5 sites avec différentes couvertures
        coverages = [10, 3, 8, 4, 6]
        for i, cov in enumerate(coverages):
            site = RNASiteVariantData(
                seqid="chr1",
                position=i * 10,
                reference=0,
                strand=1,
                coverage=cov,
                mean_quality=30.0,
                frequencies=np.array([cov-1, 1, 0, 0, 0], dtype=np.int32),
                score=0.0
            )
            counter.update(site)
        
        # Feature de 0 à 50 (50 bases)
        feature = SeqFeature(
            location=SimpleLocation(0, 50),
            id="exon1",
            type="exon"
        )
        
        total_sites = len(feature.location)
        observed_sites = counter.genome_base_freqs.sum()
        qualified_sites = counter.filtered_sites_count
        
        # Vérifications
        self.assertEqual(total_sites, 50)  # Taille de la feature
        self.assertEqual(observed_sites, 5)  # 5 observations
        self.assertEqual(qualified_sites, 3)  # 3 avec cov >= 5 (10, 8, 6)
        
        # Relations logiques
        self.assertLessEqual(observed_sites, total_sites)
        self.assertLessEqual(qualified_sites, observed_sites)


if __name__ == '__main__':
    unittest.main()
