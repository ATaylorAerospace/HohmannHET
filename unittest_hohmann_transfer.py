import pytest
import numpy as np
from matplotlib.figure import Figure
from your_module import HohmannTransfer, OrbitalConstants

@pytest.fixture(scope="session")
def constants():
    """Session-level fixture for constants to avoid repeated instantiation"""
    return OrbitalConstants()

@pytest.fixture(scope="session")
def test_cases():
    """Pre-compute common test cases to avoid repeated calculations"""
    return {
        'leo_to_meo': (300, 5000),
        'leo_to_geo': (300, 35786),
        'same_orbit': (300, 300),
        'large_transfer': (300, 100000)
    }

@pytest.fixture(scope="class")
def transfers(test_cases):
    """Create transfer objects once per test class"""
    return {name: HohmannTransfer(*altitudes) 
            for name, altitudes in test_cases.items()}

class TestHohmannTransfer:
    """Optimized test suite for HohmannTransfer class"""
    
    def test_invalid_inputs(self):
        """Batch test all invalid input cases"""
        invalid_cases = [
            (-100, 1000),
            (1000, -100),
            (0, 1000),
            (1000, 0)
        ]
        for alt1, alt2 in invalid_cases:
            with pytest.raises(ValueError):
                HohmannTransfer(alt1, alt2)

    def test_orbital_mechanics(self, transfers, constants):
        """Consolidated orbital mechanics validation"""
        for name, transfer in transfers.items():
            # Basic physical constraints
            assert transfer.r1 > constants.R_EARTH
            assert transfer.r2 > constants.R_EARTH
            assert transfer.a_transfer >= min(transfer.r1, transfer.r2)
            assert transfer.a_transfer <= max(transfer.r1, transfer.r2)
            
            # Velocity checks
            v1_calc = np.sqrt(constants.MU / transfer.r1)
            assert abs(transfer.v1 - v1_calc) < 1e-6
            
            # Energy conservation
            specific_energy = -constants.MU / (2 * transfer.a_transfer)
            assert np.isfinite(specific_energy)
            
            # Transfer time validation
            assert transfer.transfer_time > 0
            assert transfer.transfer_time < 24 * 3600  # Less than 24 hours for typical orbits

    def test_delta_v_properties(self, transfers):
        """Batch test delta-v calculations"""
        results = []
        for name, transfer in transfers.items():
            results.append({
                'name': name,
                'total_dv': transfer.total_delta_v,
                'departure_dv': transfer.delta_v_departure,
                'arrival_dv': transfer.delta_v_arrival
            })
            
            # Physical constraints
            assert transfer.total_delta_v >= 0
            assert np.isfinite(transfer.total_delta_v)
            assert transfer.total_delta_v == pytest.approx(
                transfer.delta_v_departure + transfer.delta_v_arrival
            )

    @pytest.mark.parametrize('key', ['leo_to_geo'])
    def test_known_transfers(self, transfers, key):
        """Test against known mission values"""
        transfer = transfers[key]
        if key == 'leo_to_geo':
            # Known approximate values for LEO to GEO
            assert 3.8 <= transfer.total_delta_v <= 4.0
            assert 5.0 <= transfer.transfer_time / 3600 <= 5.5  # hours

    @pytest.mark.parametrize('property', [
        'total_delta_v',
        'transfer_time',
        'a_transfer'
    ])
    def test_numerical_stability(self, transfers, property):
        """Verify numerical stability of key properties"""
        for transfer in transfers.values():
            value = getattr(transfer, property)
            assert np.isfinite(value)
            assert not np.isnan(value)
            assert not np.isinf(value)

    def test_visualization(self, transfers):
        """Optimized visualization testing"""
        transfer = next(iter(transfers.values()))  # Test only one case
        fig = transfer.visualize_transfer()
        
        # Basic figure validation
        assert isinstance(fig, Figure)
        ax = fig.axes[0]
        
        # Efficient plot element validation
        children = ax.get_children()
        has_surface = any(isinstance(child, plt.Surface) for child in children)
        has_lines = len([c for c in children if isinstance(c, plt.Line2D)]) >= 2
        
        assert has_surface  # Earth representation
        assert has_lines    # Orbit lines
        assert ax.get_xlabel() and ax.get_ylabel() and ax.get_zlabel()

    def test_same_orbit_case(self, transfers):
        """Special case: same orbit transfer"""
        transfer = transfers['same_orbit']
        assert transfer.total_delta_v < 1e-6
        assert transfer.delta_v_departure < 1e-6
        assert transfer.delta_v_arrival < 1e-6

    def test_large_transfer_stability(self, transfers):
        """Test numerical stability for large transfers"""
        transfer = transfers['large_transfer']
        assert np.isfinite(transfer.total_delta_v)
        assert np.isfinite(transfer.transfer_time)
        assert transfer.total_delta_v > 0
