# import numpy as np
class Base:
    def __init__(self, base_id, base_type, prev_id, next_id, position,
                 backbone, normal, velocity, angular_velocity, strand_id = -1):
        self.base_id = base_id
        self.base_type = base_type
        self.prev_id = prev_id
        self.next_id = next_id
        self.position = position
        self.backbone = backbone
        self.normal = normal
        self.velocity = velocity
        self.angular_velocity = angular_velocity
        self.strand_id = strand_id
    
    def set_position(self, position):
        assert len(position) == 3
        self.position = position
        return self.position

    @staticmethod
    def parse_string(s):
        pass

    @staticmethod
    def parse_list(params):
        (
            position,
            backbone,
            normal,
            velocity,
            angular_velocity,
            base_id,
            strand_id,
            base_type,
            prev_id,
            next_id,
        ) = params
        return Base(base_id, base_type, prev_id, next_id, position,
                    backbone, normal, velocity, angular_velocity, strand_id)
