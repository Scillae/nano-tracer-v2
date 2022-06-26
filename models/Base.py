# import numpy as np
class Base:
    def __init__(self, base_id, base_type, prev_id, next_id, position,
                 backbone, normal, velocity, angular_velocity, strand_id=-1):
        """
        init Base
        :param base_id:
        :param base_type:
        :param prev_id:
        :param next_id:
        :param position:
        :param backbone:
        :param normal:
        :param velocity:
        :param angular_velocity:
        :param strand_id:
        """
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
        """
        set position of base
        :param position:
        :return:
        """
        assert len(position) == 3
        self.position = position
        return self.position

    @staticmethod
    def parse_string(s):
        pass

    @staticmethod
    def parse_list(params):
        """
        unpack list/tuple to Base instance
        :param params: list/tuple of parameters
        :return: Base instance
        """
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
